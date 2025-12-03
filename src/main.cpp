#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>
#include <map>
#include <deque>
#include <dirent.h>

// --- Constants ---
int SCREEN_WIDTH = 1200;
int SCREEN_HEIGHT = 800;
const float ATOM_RADIUS = 12.0f;
const float MOLECULE_SCALE = 45.0f;
const float MAX_SPEED = 2.0f;

// Bond lengths in Angstroms (from PubChem/CRC Handbook)
// Converted to pixels: 1 Angstrom = ~40px for visual clarity
const float BOND_SINGLE = 54.0f;   // C-C ~1.54A
const float BOND_DOUBLE = 49.2f;  // C=C ~1.34A (shorter)
const float BOND_TRIPLE = 44.0f;  // C#C ~1.20A (even shorter)
const float BOND_AROMATIC = 51.6f; // Benzene ~1.39A
const float BOND_HYDROGEN_LEN = 70.0f; // ~1.8-2.0A

// --- Enums & Structs ---
enum AppState { STATE_LOADING, STATE_SIMULATION };
enum ComplexityLevel { LEVEL_SKELETAL = 1, LEVEL_ORGANIC = 2, LEVEL_LEWIS = 3 };

struct Element {
    std::string symbol;
    int valency;
    int valence_electrons; // For Lewis structures
    SDL_Color color;
};

struct Bond {
    int target_id;
    int order;
};

struct Atom {
    float x, y;
    float vx, vy;
    std::string symbol;
    std::string parent_molecule; // For specificity (e.g. "Adenine")
    int max_bonds;
    int valence_electrons;
    int current_valency; // Actual valency used (Single=1, Double=2)
    SDL_Color color;
    int id;
    std::vector<Bond> bonds;
};

struct MoleculeTemplate {
    std::string name;
    struct TAtom { int id; float x, y; std::string symbol; };
    struct TBond { int id1, id2; int type; };
    std::vector<TAtom> t_atoms;
    std::vector<TBond> t_bonds;
};

// --- Globals ---
std::vector<Atom> atoms;
std::vector<MoleculeTemplate> templates;
std::deque<std::string> system_logs;
int global_atom_id = 0;
int current_template_index = 0;
AppState current_state = STATE_LOADING;
ComplexityLevel complexity_level = LEVEL_LEWIS;
bool is_fullscreen = false;
float loading_progress = 0.0f;
int last_mouse_x = 0, last_mouse_y = 0;

// --- Data ---
// Covalent radii in pm (picometers) - used for spring target distance calculation
// Source: PubChem Periodic Table
struct AtomData {
    std::string symbol;
    int valency;
    int valence_electrons;
    float covalent_radius; // in pm
    SDL_Color color;
};

AtomData getAtomData(std::string sym) {
    // Symbol, Valency, ValenceE, CovalentRadius(pm), Color (CPK-ish)
    if (sym == "C") return {"C", 4, 4, 77.0f, {50, 50, 50, 255}};       // Carbon: Grey
    if (sym == "O") return {"O", 2, 6, 73.0f, {200, 50, 50, 255}};      // Oxygen: Red
    if (sym == "N") return {"N", 3, 5, 75.0f, {50, 50, 200, 255}};      // Nitrogen: Blue
    if (sym == "H") return {"H", 1, 1, 32.0f, {200, 200, 200, 255}};    // Hydrogen: White
    if (sym == "S") return {"S", 2, 6, 102.0f, {200, 200, 50, 255}};    // Sulfur: Yellow
    if (sym == "P") return {"P", 3, 5, 106.0f, {255, 128, 0, 255}};     // Phosphorus: Orange
    if (sym == "Cl") return {"Cl", 1, 7, 99.0f, {50, 200, 50, 255}};    // Chlorine: Green
    if (sym == "Br") return {"Br", 1, 7, 114.0f, {150, 50, 50, 255}};   // Bromine: Dark Red
    if (sym == "F") return {"F", 1, 7, 64.0f, {144, 224, 80, 255}};     // Fluorine: Light Green
    return {"?", 0, 0, 70.0f, {100, 0, 100, 255}};
}

// Legacy wrapper for compatibility
Element getElement(std::string sym) {
    AtomData d = getAtomData(sym);
    return {d.symbol, d.valency, d.valence_electrons, d.color};
}

// Calculate ideal bond length between two atoms (sum of covalent radii, scaled)
float getBondLength(std::string sym1, std::string sym2, int bondOrder = 1) {
    AtomData a1 = getAtomData(sym1);
    AtomData a2 = getAtomData(sym2);
    // Sum of covalent radii in pm, scaled to pixels (0.35 px/pm)
    float baseLen = (a1.covalent_radius + a2.covalent_radius) * 0.35f;
    // Shorten for multiple bonds
    if (bondOrder == 2) baseLen *= 0.87f; // ~13% shorter
    if (bondOrder == 3) baseLen *= 0.78f; // ~22% shorter
    if (bondOrder == 4) return BOND_HYDROGEN_LEN; // Hydrogen Bond
    return baseLen;
}

// --- Helper Functions ---
void log(std::string msg) {
    system_logs.push_back(msg);
    if (system_logs.size() > 100) system_logs.pop_front(); // Increased log buffer to 100
}

float distSq(Atom& a, Atom& b) {
    return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
}

void spawnAtom(float x, float y, std::string sym, std::string parent = "") {
    Element e = getElement(sym);
    Atom a;
    a.x = x; a.y = y;
    a.vx = ((float)(rand() % 100) / 50.0f - 1.0f) * MAX_SPEED;
    a.vy = ((float)(rand() % 100) / 50.0f - 1.0f) * MAX_SPEED;
    a.symbol = e.symbol;
    a.parent_molecule = parent;
    a.max_bonds = e.valency;
    a.valence_electrons = e.valence_electrons;
    a.current_valency = 0;
    a.color = e.color;
    a.id = global_atom_id++;
    atoms.push_back(a);
}

void loadMolecule(std::string filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        log("[ERROR] Failed to open " + filepath);
        return;
    }

    MoleculeTemplate templ;
    std::string line;
    log("[SYSTEM] Accessing " + filepath + "...");
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            if (line.find("Molecule:") != std::string::npos) {
                templ.name = line.substr(line.find(":") + 2);
                log("[FOUND] Compound: " + templ.name);
            }
            continue;
        }
        std::stringstream ss(line);
        std::string type;
        ss >> type;
        
        if (type == "ATOM") {
            int id; float x, y; std::string sym;
            ss >> id >> x >> y >> sym;
            templ.t_atoms.push_back({id, x, y, sym});
        }
        else if (type == "BOND") {
            int id1, id2, btype;
            ss >> id1 >> id2 >> btype;
            templ.t_bonds.push_back({id1, id2, btype});
        }
    }
    templates.push_back(templ);
    log("[SUCCESS] Compiled " + templ.name);
}

void loadAllMolecules(std::string dirPath) {
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(dirPath.c_str())) != NULL) {
        std::vector<std::string> files;
        while ((ent = readdir(dir)) != NULL) {
            std::string fname = ent->d_name;
            if (fname.length() > 5 && fname.substr(fname.length()-5) == ".chem") {
                files.push_back(dirPath + "/" + fname);
            }
        }
        closedir(dir);
        
        // Sort to ensure deterministic order
        for(const auto& f : files) {
            loadMolecule(f);
        }
    } else {
        log("[ERROR] Could not open " + dirPath);
    }
}

void spawnMolecule(int templateIdx, float cx, float cy) {
    if (templateIdx < 0 || templateIdx >= (int)templates.size()) return;
    MoleculeTemplate& t = templates[templateIdx];
    std::map<int, int> id_map;
    
    for (auto& ta : t.t_atoms) {
        float x = cx + ta.x * MOLECULE_SCALE;
        float y = cy - ta.y * MOLECULE_SCALE;
        spawnAtom(x, y, ta.symbol, t.name);
        id_map[ta.id] = atoms.back().id;
        atoms.back().vx = 0; atoms.back().vy = 0;
    }
    
    for (auto& tb : t.t_bonds) {
        int gid1 = id_map[tb.id1];
        int gid2 = id_map[tb.id2];
        Atom* a1 = nullptr; Atom* a2 = nullptr;
        for (auto& a : atoms) {
            if (a.id == gid1) a1 = &a;
            if (a.id == gid2) a2 = &a;
        }
        if (a1 && a2) {
            a1->bonds.push_back({a2->id, tb.type});
            a1->current_valency += tb.type;
            a2->bonds.push_back({a1->id, tb.type});
            a2->current_valency += tb.type;
        }
    }
}

// --- Rendering Helpers ---

void drawText(SDL_Renderer* r, TTF_Font* f, std::string text, int x, int y, SDL_Color c) {
    if (!f) return;
    SDL_Surface* surf = TTF_RenderText_Solid(f, text.c_str(), c);
    if (surf) {
        SDL_Texture* tex = SDL_CreateTextureFromSurface(r, surf);
        SDL_Rect rect = {x, y, surf->w, surf->h};
        SDL_RenderCopy(r, tex, NULL, &rect);
        SDL_FreeSurface(surf);
        SDL_DestroyTexture(tex);
    }
}

void drawLonePairs(SDL_Renderer* r, Atom& a) {
    // Calculate lone pairs: (Valence - Bonds) / 2
    // Note: This is a simplification. Charge is assumed 0.
    int lone_pairs = (a.valence_electrons - a.current_valency) / 2;
    if (lone_pairs <= 0) return;

    // Calculate vector sum of bonds to find "empty" side
    float sumX = 0, sumY = 0;
    for (auto& b : a.bonds) {
        for (auto& other : atoms) {
            if (other.id == b.target_id) {
                float dx = other.x - a.x;
                float dy = other.y - a.y;
                float len = sqrt(dx*dx + dy*dy);
                if (len > 0) { sumX += dx/len; sumY += dy/len; }
            }
        }
    }

    // Normalize repulsion vector
    float repLen = sqrt(sumX*sumX + sumY*sumY);
    float dirX = (repLen > 0) ? -(sumX/repLen) : 1.0f; // Opposite to bonds
    float dirY = (repLen > 0) ? -(sumY/repLen) : 0.0f;

    // Draw pairs
    SDL_SetRenderDrawColor(r, 255, 255, 0, 200); // Yellow electrons
    float dist = ATOM_RADIUS + 5.0f;
    
    if (lone_pairs == 1) {
        SDL_RenderDrawPoint(r, a.x + dirX*dist - 2, a.y + dirY*dist);
        SDL_RenderDrawPoint(r, a.x + dirX*dist + 2, a.y + dirY*dist);
    } else if (lone_pairs == 2) {
        // Split angle
        float angle = atan2(dirY, dirX);
        float a1 = angle - 0.5f;
        float a2 = angle + 0.5f;
        
        SDL_RenderDrawPoint(r, a.x + cos(a1)*dist, a.y + sin(a1)*dist);
        SDL_RenderDrawPoint(r, a.x + cos(a1)*dist+2, a.y + sin(a1)*dist);
        
        SDL_RenderDrawPoint(r, a.x + cos(a2)*dist, a.y + sin(a2)*dist);
        SDL_RenderDrawPoint(r, a.x + cos(a2)*dist+2, a.y + sin(a2)*dist);
    }
}

// --- Main ---

int main(int argc, char* args[]) {
    srand(time(NULL));

    if (SDL_Init(SDL_INIT_VIDEO) < 0) return 1;
    if (TTF_Init() == -1) return 1;

    // Borderless window for custom controls
    SDL_Window* window = SDL_CreateWindow("ChemGen", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
                                        SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN | SDL_WINDOW_BORDERLESS);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    TTF_Font* font = TTF_OpenFont("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 14);
    if (!font) font = TTF_OpenFont("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 14);
    TTF_Font* fontLarge = TTF_OpenFont("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 24);
    if (!fontLarge) fontLarge = TTF_OpenFont("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 24);

    bool quit = false;
    SDL_Event e;
    
    // Loading Sequence Data
    log("[INIT] Molecular Art Station v1.0");
    log("[MEM] Allocating particle buffers...");
    
    // Simulation Loop
    while (!quit) {
        // Update window size for relative UI
        int winW, winH;
        SDL_GetWindowSize(window, &winW, &winH);
        SCREEN_WIDTH = winW;
        SCREEN_HEIGHT = winH;

        int mouseX, mouseY;
        Uint32 mouseState = SDL_GetMouseState(&mouseX, &mouseY);

        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) quit = true;
            if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_ESCAPE) quit = true;
                if (e.key.keysym.sym == SDLK_SPACE) atoms.clear();
                if (e.key.keysym.sym == SDLK_r) {
                    atoms.clear();
                    log("[SYS] Reset Simulation");
                }
                if (e.key.keysym.sym == SDLK_f) {
                    is_fullscreen = !is_fullscreen;
                    SDL_SetWindowFullscreen(window, is_fullscreen ? SDL_WINDOW_FULLSCREEN_DESKTOP : 0);
                }
                // Number keys for complexity
                if (e.key.keysym.sym == SDLK_1 || e.key.keysym.sym == SDLK_KP_1) complexity_level = LEVEL_SKELETAL;
                if (e.key.keysym.sym == SDLK_2 || e.key.keysym.sym == SDLK_KP_2) complexity_level = LEVEL_ORGANIC;
                if (e.key.keysym.sym == SDLK_3 || e.key.keysym.sym == SDLK_KP_3) complexity_level = LEVEL_LEWIS;
                
                // Enter for Spawn & Iterate
                if (e.key.keysym.sym == SDLK_RETURN || e.key.keysym.sym == SDLK_RETURN2 || e.key.keysym.sym == SDLK_KP_ENTER) {
                    spawnMolecule(current_template_index, (float)mouseX, (float)mouseY);
                    // Iterate to next template
                    if (!templates.empty()) {
                        current_template_index = (current_template_index + 1) % templates.size();
                        log("[SEL] " + templates[current_template_index].name);
                    }
                }
            }
            // UI Clicks
            if (e.type == SDL_MOUSEBUTTONDOWN) {
                // Window Controls (Top Right)
                if (mouseX > SCREEN_WIDTH - 30 && mouseY < 30) quit = true; // Close
                if (mouseX > SCREEN_WIDTH - 60 && mouseX < SCREEN_WIDTH - 30 && mouseY < 30) { // Maximize/Restore
                    is_fullscreen = !is_fullscreen;
                    SDL_SetWindowFullscreen(window, is_fullscreen ? SDL_WINDOW_FULLSCREEN_DESKTOP : 0);
                }
                if (mouseX > SCREEN_WIDTH - 90 && mouseX < SCREEN_WIDTH - 60 && mouseY < 30) { // Minimize
                    SDL_MinimizeWindow(window);
                }
                
                // Complexity Slider (Bottom Left)
                if (mouseY > SCREEN_HEIGHT - 40 && mouseY < SCREEN_HEIGHT - 10) {
                    if (mouseX > 10 && mouseX < 110) complexity_level = LEVEL_SKELETAL;
                    if (mouseX > 120 && mouseX < 220) complexity_level = LEVEL_ORGANIC;
                    if (mouseX > 230 && mouseX < 330) complexity_level = LEVEL_LEWIS;
                    if (mouseX > 330 && mouseX < 410) { // Reset
                        atoms.clear();
                        log("[SYS] Reset Simulation");
                    }
                    // Mute spawn on click
                    continue; 
                }
            }
        }

        // --- LOGIC ---
        
        if (current_state == STATE_LOADING) {
            loading_progress += 0.5f;
            if (loading_progress == 20.0f) {
                // Load all .chem files from data/
                loadAllMolecules("data");
                if (templates.empty()) {
                    // Fallback if no files found
                    log("[WARN] No .chem files found in data/");
                }
            }
            if (loading_progress > 60.0f) {
                current_state = STATE_SIMULATION;
                log("[SYS] Simulation Started.");
            }
        } else {
            // Simulation Logic
            if (mouseState & SDL_BUTTON(SDL_BUTTON_LEFT) && mouseY > 40 && mouseY < SCREEN_HEIGHT - 50) { // Don't spawn on UI areas
                int r = rand() % 100;
                std::string sym = "C";
                // Adjusted probabilities: Less H, more heteroatoms
                if (r > 50) sym = "O"; 
                if (r > 75) sym = "N";
                // Add S and P with lower probability
                if (r > 90) sym = "S"; 
                if (r > 95) sym = "P";
                // Very rare H spawn (mostly H should come from filling valency, not random spawn)
                if (r > 98) sym = "H"; 

                spawnAtom((float)mouseX + (rand()%20-10), (float)mouseY + (rand()%20-10), sym);
            }
            if (mouseState & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
                if (rand() % 10 == 0) {
                    spawnMolecule(current_template_index, (float)mouseX, (float)mouseY);
                    current_template_index = (current_template_index + 1) % templates.size();
                }
            }

            // --- PHYSICS ENGINE REWRITE ---
            
            // 1. Integration & Wall Constraints
            for (auto& a : atoms) {
                a.x += a.vx; 
                a.y += a.vy;
                
                // Hard Wall Bounce
                if (a.x < ATOM_RADIUS) { a.x = ATOM_RADIUS; a.vx = abs(a.vx) * 0.5f; }
                if (a.x > SCREEN_WIDTH - ATOM_RADIUS) { a.x = SCREEN_WIDTH - ATOM_RADIUS; a.vx = -abs(a.vx) * 0.5f; }
                if (a.y < ATOM_RADIUS) { a.y = ATOM_RADIUS; a.vy = abs(a.vy) * 0.5f; }
                if (a.y > SCREEN_HEIGHT - ATOM_RADIUS) { a.y = SCREEN_HEIGHT - ATOM_RADIUS; a.vy = -abs(a.vy) * 0.5f; }
                
                // Friction
                a.vx *= 0.98f; 
                a.vy *= 0.98f;
            }

            // 2. Dynamic Bonding
            for (size_t i = 0; i < atoms.size(); i++) {
                if (atoms[i].current_valency >= atoms[i].max_bonds) continue;
                
                for (size_t j = i + 1; j < atoms.size(); j++) {
                    if (atoms[j].current_valency >= atoms[j].max_bonds) continue;

                    // Check if already bonded
                    bool already_bonded = false;
                    for (auto& b : atoms[i].bonds) if (b.target_id == atoms[j].id) already_bonded = true;
                    if (already_bonded) continue;

                    float d2 = distSq(atoms[i], atoms[j]);
                    // Snap distance: slightly larger than single bond to catch pairs
                    float snapDist = BOND_SINGLE * 1.3f; 
                    
                    if (d2 < snapDist*snapDist) {
                        // Logic: Only bond if compatible
                        // Prevent random H-H or H-C unless specific conditions (optional)
                        // Prevent random H spawning
                        
                        atoms[i].bonds.push_back({atoms[j].id, 1});
                        atoms[i].current_valency++;
                        atoms[j].bonds.push_back({atoms[i].id, 1});
                        atoms[j].current_valency++;
                    }
                }
            }

            // 2b. Hydrogen Bonding (Intermolecular)
            // Allows "saturated" atoms to bond if they are H-bond donors/acceptors
            // Broadened rules to allow non-canonical pairing (Hoogsteen, etc.) and backbone interactions
            for (size_t i = 0; i < atoms.size(); i++) {
                for (size_t j = i + 1; j < atoms.size(); j++) {
                    // Check for H-bond candidates: H (on N/O/F/P) <-> N/O/F/P
                    // We treat N, O, F, P as electronegative atoms capable of H-bonding
                    bool isHBond = false;
                    
                    std::string s1 = atoms[i].symbol;
                    std::string s2 = atoms[j].symbol;
                    
                    bool a1_is_H = (s1 == "H");
                    bool a2_is_H = (s2 == "H");
                    
                    if ((a1_is_H || a2_is_H) && !(a1_is_H && a2_is_H)) {
                        Atom* hAtom = a1_is_H ? &atoms[i] : &atoms[j];
                        Atom* accAtom = a1_is_H ? &atoms[j] : &atoms[i];
                        
                        // Check Acceptor (N, O, F, P, Cl)
                        std::string accSym = accAtom->symbol;
                        if (accSym == "N" || accSym == "O" || accSym == "F" || accSym == "P" || accSym == "Cl") {
                            
                            // Check if H is polar (bonded to N, O, F, P, Cl)
                            bool hIsPolar = false;
                            for(auto& b : hAtom->bonds) {
                                for(auto& other : atoms) {
                                    if(other.id == b.target_id) {
                                        std::string dSym = other.symbol;
                                        if(dSym == "N" || dSym == "O" || dSym == "F" || dSym == "P" || dSym == "Cl") {
                                            hIsPolar = true;
                                        }
                                        break;
                                    }
                                }
                            }
                            
                            if (hIsPolar) {
                                // Check distance for H-bond
                                float d2 = distSq(atoms[i], atoms[j]);
                                // Allow bonding if close enough
                                if (d2 < (BOND_HYDROGEN_LEN * 1.2f)*(BOND_HYDROGEN_LEN * 1.2f)) {
                                    // Check if already bonded
                                    bool already_bonded = false;
                                    for(auto& b : atoms[i].bonds) if(b.target_id == atoms[j].id) already_bonded = true;
                                    
                                    if (!already_bonded) {
                                        // Specificity Check: Only allow Canonical Base Pairs
                                        // A-T (or A-U), G-C
                                        std::string p1 = atoms[i].parent_molecule;
                                        std::string p2 = atoms[j].parent_molecule;
                                        
                                        // If parents are empty (random atoms), allow generic H-bonds
                                        // If parents are known nucleotides, enforce rules
                                        bool allowed = true;
                                        if (!p1.empty() && !p2.empty()) {
                                            // Normalize names (remove extension if present, though loadMolecule strips it usually)
                                            // Check pairs
                                            bool pairAT = (p1 == "Adenine" && (p2 == "Thymine" || p2 == "Uracil")) || 
                                                          (p2 == "Adenine" && (p1 == "Thymine" || p1 == "Uracil"));
                                            bool pairGC = (p1 == "Guanine" && p2 == "Cytosine") || 
                                                          (p2 == "Guanine" && p1 == "Cytosine");
                                            
                                            if (!pairAT && !pairGC) allowed = false;
                                        }

                                        if (allowed) isHBond = true;
                                    }
                                }
                            }
                        }
                    }

                    if (isHBond) {
                         atoms[i].bonds.push_back({atoms[j].id, 4});
                         atoms[j].bonds.push_back({atoms[i].id, 4});
                    }
                }
            }

            // 3. Hard Lock Constraints (Iterative Solver)
            // This replaces springs with rigid distance constraints
            int iterations = 10; // More iterations = stiffer bonds
            for (int k = 0; k < iterations; k++) {
                for (size_t i = 0; i < atoms.size(); i++) {
                    for (auto& b : atoms[i].bonds) {
                        int tid = b.target_id;
                        // Find neighbor
                        int nIdx = -1;
                        for(size_t m=0; m<atoms.size(); m++) if(atoms[m].id == tid) { nIdx = m; break; }
                        if (nIdx == -1) continue;

                        // Avoid double processing
                        if (atoms[i].id > atoms[nIdx].id) continue;

                        Atom& a1 = atoms[i];
                        Atom& a2 = atoms[nIdx];

                        float dx = a2.x - a1.x;
                        float dy = a2.y - a1.y;
                        float dist = sqrt(dx*dx + dy*dy);
                        if (dist < 0.001f) dist = 0.001f;

                        // Determine Bond Target Length
                        // Use stored bond order!
                        float target = getBondLength(a1.symbol, a2.symbol, b.order);

                        // Constraint Correction
                        float diff = dist - target;
                        float correction = diff / dist * 0.5f; // Split move equally
                        
                        float offX = dx * correction;
                        float offY = dy * correction;

                        a1.x += offX; a1.y += offY;
                        a2.x -= offX; a2.y -= offY;
                        
                        // Update velocities to reflect the constraint (optional but helps stability)
                        // For now, position correction is enough for "Hard Lock"
                    }
                }
            }

            // 4. Molecule Collision (Bouncing entire molecules)
            // Step A: Identify Molecules (Connected Components)
            std::vector<std::vector<int>> molecules;
            std::vector<bool> visited(atoms.size(), false);
            std::map<int, int> idToIndex;
            for(size_t i=0; i<atoms.size(); i++) idToIndex[atoms[i].id] = i;

            for(size_t i=0; i<atoms.size(); i++) {
                if(visited[i]) continue;
                std::vector<int> mol;
                std::deque<int> q;
                q.push_back(i);
                visited[i] = true;
                while(!q.empty()) {
                    int curr = q.front(); q.pop_front();
                    mol.push_back(curr);
                    for(auto& b : atoms[curr].bonds) {
                        if(idToIndex.count(b.target_id)) {
                            int neighborIdx = idToIndex[b.target_id];
                            if(!visited[neighborIdx]) {
                                visited[neighborIdx] = true;
                                q.push_back(neighborIdx);
                            }
                        }
                    }
                }
                molecules.push_back(mol);
            }

            // Step B: Collide Molecules
            // REPLACED Bounding Circle with Atom-Atom Repulsion for better shape fitting
            // This allows irregular molecules (like nucleotides) to interlock and bond
            for (size_t i = 0; i < molecules.size(); i++) {
                for (size_t j = i + 1; j < molecules.size(); j++) {
                    // Detailed Atom-Atom check between molecules
                    for(int idx1 : molecules[i]) {
                        for(int idx2 : molecules[j]) {
                            float dx = atoms[idx1].x - atoms[idx2].x;
                            float dy = atoms[idx1].y - atoms[idx2].y;
                            float d2 = dx*dx + dy*dy;
                            float minR = ATOM_RADIUS * 2.0f; // Simple hard shell
                            
                            if (d2 < minR*minR && d2 > 0) {
                                float dist = sqrt(d2);
                                float overlap = minR - dist;
                                float nx = dx/dist;
                                float ny = dy/dist;
                                
                                // Push apart
                                atoms[idx1].x += nx * overlap * 0.5f;
                                atoms[idx1].y += ny * overlap * 0.5f;
                                atoms[idx2].x -= nx * overlap * 0.5f;
                                atoms[idx2].y -= ny * overlap * 0.5f;
                                
                                // Transfer momentum (elastic-ish)
                                float v1n = atoms[idx1].vx * nx + atoms[idx1].vy * ny;
                                float v2n = atoms[idx2].vx * nx + atoms[idx2].vy * ny;
                                
                                // Simple exchange with damping
                                atoms[idx1].vx -= (v1n - v2n) * nx * 0.8f;
                                atoms[idx1].vy -= (v1n - v2n) * ny * 0.8f;
                                atoms[idx2].vx += (v1n - v2n) * nx * 0.8f;
                                atoms[idx2].vy += (v1n - v2n) * ny * 0.8f;
                            }
                        }
                    }
                }
            }
            // 5. Auto-Saturation (Lewis Mode)
            // If in Lewis mode, ensure all atoms have full valency by spawning Hydrogens
            if (complexity_level == LEVEL_LEWIS) {
                static int saturation_timer = 0;
                saturation_timer++;
                if (saturation_timer > 30) { // Check every ~0.5s
                    saturation_timer = 0;
                    // Iterate backwards to allow spawning without invalidating iterators immediately
                    // (Though vector reallocation might still happen, so we use indices)
                    size_t current_size = atoms.size(); 
                    for (size_t i = 0; i < current_size; i++) {
                        // Only saturate if the atom is somewhat stable (not flying too fast)
                        if (abs(atoms[i].vx) > 0.5f || abs(atoms[i].vy) > 0.5f) continue;

                        int missing = atoms[i].max_bonds - atoms[i].current_valency;
                        if (missing > 0) {
                            // Calculate base angle based on existing bonds to avoid overlap
                            float base_angle = 0.0f;
                            if (!atoms[i].bonds.empty()) {
                                // Find angle away from first bond to distribute H on the other side
                                for(auto& other : atoms) {
                                    if(other.id == atoms[i].bonds[0].target_id) {
                                        base_angle = atan2(atoms[i].y - other.y, atoms[i].x - other.x);
                                        break;
                                    }
                                }
                            } else {
                                base_angle = (float)(rand() % 360) * 3.14159f / 180.0f;
                            }

                            // Spawn all needed Hydrogens at once
                            for (int k = 0; k < missing; k++) {
                                // Spread hydrogens out in a fan or circle
                                float angle_offset = (missing == 1) ? 0.0f : ((k + 1) * (3.14159f / (missing + 1)) - 1.57f);
                                if (atoms[i].bonds.empty()) angle_offset = k * (6.28f / missing); // Full circle if no bonds
                                
                                float angle = base_angle + angle_offset;
                                // Add slight randomness for organic feel
                                angle += ((rand()%100)/100.0f - 0.5f) * 0.2f;

                                float dist = BOND_SINGLE;
                                float hx = atoms[i].x + cos(angle) * dist;
                                float hy = atoms[i].y + sin(angle) * dist;
                                
                                spawnAtom(hx, hy, "H");
                                // Force bond immediately
                                Atom& h = atoms.back();
                                h.vx = atoms[i].vx; h.vy = atoms[i].vy; // Match velocity
                                
                                atoms[i].bonds.push_back({h.id, 1});
                                atoms[i].current_valency++;
                                h.bonds.push_back({atoms[i].id, 1});
                                h.current_valency++;
                            }
                        }
                    }
                }
            }
        }

        // --- RENDER ---
        SDL_SetRenderDrawColor(renderer, 20, 20, 25, 255);
        SDL_RenderClear(renderer);

        if (current_state == STATE_LOADING) {
            // Draw Logs
            int y = 100;
            for (auto& log : system_logs) {
                drawText(renderer, font, log, 50, y, {0, 255, 0, 255});
                y += 20;
            }
            // Draw Loading Bar
            SDL_Rect bar = {50, y + 20, (int)(loading_progress * 5), 10};
            SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
            SDL_RenderFillRect(renderer, &bar);
            
            // Draw Preview Molecule (Rotating)
            if (!templates.empty()) {
                MoleculeTemplate& t = templates.back();
                float cx = SCREEN_WIDTH / 2 + 200;
                float cy = SCREEN_HEIGHT / 2;
                float angle = loading_progress * 0.1f;
                
                for (auto& ta : t.t_atoms) {
                    // Rotate
                    float rx = ta.x * cos(angle) - ta.y * sin(angle);
                    float ry = ta.x * sin(angle) + ta.y * cos(angle);
                    
                    SDL_Rect r = {(int)(cx + rx*40), (int)(cy + ry*40), 10, 10};
                    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
                    SDL_RenderFillRect(renderer, &r);
                }
            }
        } 
        else {
            // Draw Bonds
            SDL_SetRenderDrawColor(renderer, 100, 100, 100, 255);
            for (auto& a : atoms) {
                for (auto& b : a.bonds) {
                    int tid = b.target_id;
                    // Avoid drawing twice
                    if (a.id > tid) continue;

                    for (auto& t : atoms) {
                        if (t.id == tid) {
                            int bondOrder = b.order;

                            // Draw Lines
                            if (bondOrder == 2) {
                                SDL_RenderDrawLine(renderer, a.x-3, a.y-3, t.x-3, t.y-3);
                                SDL_RenderDrawLine(renderer, a.x+3, a.y+3, t.x+3, t.y+3);
                            } else if (bondOrder == 3) {
                                SDL_RenderDrawLine(renderer, a.x, a.y, t.x, t.y);
                                SDL_RenderDrawLine(renderer, a.x-4, a.y-4, t.x-4, t.y-4);
                                SDL_RenderDrawLine(renderer, a.x+4, a.y+4, t.x+4, t.y+4);
                            } else if (bondOrder == 4) {
                                // Hydrogen Bond (Dashed Line - simulated with dots)
                                float dx = t.x - a.x;
                                float dy = t.y - a.y;
                                float dist = sqrt(dx*dx + dy*dy);
                                int segments = (int)(dist / 5.0f);
                                SDL_SetRenderDrawColor(renderer, 100, 100, 255, 150); // Light Blue
                                for(int k=0; k<segments; k+=2) {
                                    float x1 = a.x + (dx/dist) * (k*5.0f);
                                    float y1 = a.y + (dy/dist) * (k*5.0f);
                                    float x2 = a.x + (dx/dist) * ((k+1)*5.0f);
                                    float y2 = a.y + (dy/dist) * ((k+1)*5.0f);
                                    SDL_RenderDrawLine(renderer, (int)x1, (int)y1, (int)x2, (int)y2);
                                }
                                SDL_SetRenderDrawColor(renderer, 100, 100, 100, 255); // Reset
                            } else {
                                SDL_RenderDrawLine(renderer, a.x, a.y, t.x, t.y);
                            }
                            break;
                        }
                    }
                }
            }

            // Draw Atoms
            for (auto& a : atoms) {
                bool draw = true;
                bool drawTextSym = true;

                if (complexity_level == LEVEL_SKELETAL) {
                    draw = false; // Hide atoms in skeletal
                }
                else if (complexity_level == LEVEL_ORGANIC) {
                    // Hide Hydrogens on Carbons
                    if (a.symbol == "H") {
                        for (auto& b : a.bonds) {
                            for (auto& n : atoms) if (n.id == b.target_id && n.symbol == "C") draw = false;
                        }
                    }
                    // Carbons are just dots
                    if (a.symbol == "C") drawTextSym = false;
                }

                if (draw) {
                    SDL_Rect rect = {(int)(a.x - ATOM_RADIUS), (int)(a.y - ATOM_RADIUS), (int)(ATOM_RADIUS*2), (int)(ATOM_RADIUS*2)};
                    SDL_SetRenderDrawColor(renderer, a.color.r, a.color.g, a.color.b, a.color.a);
                    if (complexity_level == LEVEL_ORGANIC && a.symbol == "C") {
                        // Small dot for Carbon
                        SDL_Rect small = {(int)a.x-3, (int)a.y-3, 6, 6};
                        SDL_RenderFillRect(renderer, &small);
                    } else {
                        SDL_RenderFillRect(renderer, &rect);
                        if (drawTextSym) {
                            SDL_Color tc = {255, 255, 255};
                            if (a.symbol == "H") tc = {0,0,0};
                            drawText(renderer, font, a.symbol, a.x-5, a.y-8, tc);
                        }
                    }
                }
                
                if (complexity_level == LEVEL_LEWIS) {
                    drawLonePairs(renderer, a);
                }
            }
            
            // Draw UI Overlay
            // Get current window size for relative positioning
            int winW, winH;
            SDL_GetWindowSize(window, &winW, &winH);
            SCREEN_WIDTH = winW;
            SCREEN_HEIGHT = winH;

            // Window Controls (Top Right) - Line icons only, appear on hover
            if (mouseY < 35) {
                SDL_SetRenderDrawColor(renderer, 30, 30, 35, 180);
                SDL_Rect header = {winW - 70, 0, 70, 30};
                SDL_RenderFillRect(renderer, &header);
                
                // Maximize/Restore icon (two overlapping squares or single square)
                SDL_SetRenderDrawColor(renderer, 180, 180, 180, 255);
                if (is_fullscreen) {
                    // Restore icon: two overlapping rectangles
                    SDL_RenderDrawLine(renderer, winW-55, 8, winW-45, 8);
                    SDL_RenderDrawLine(renderer, winW-55, 8, winW-55, 18);
                    SDL_RenderDrawLine(renderer, winW-58, 12, winW-48, 12);
                    SDL_RenderDrawLine(renderer, winW-58, 12, winW-58, 22);
                    SDL_RenderDrawLine(renderer, winW-58, 22, winW-48, 22);
                    SDL_RenderDrawLine(renderer, winW-48, 12, winW-48, 22);
                } else {
                    // Maximize icon: single rectangle
                    SDL_Rect maxIcon = {winW-58, 8, 14, 14};
                    SDL_RenderDrawRect(renderer, &maxIcon);
                }
                
                // Close icon (X)
                SDL_SetRenderDrawColor(renderer, 255, 80, 80, 255);
                SDL_RenderDrawLine(renderer, winW-25, 8, winW-15, 22);
                SDL_RenderDrawLine(renderer, winW-25, 22, winW-15, 8);
            }
            
            // Complexity Slider (Bottom Left) - Relative to window
            int by = winH - 35;
            SDL_SetRenderDrawColor(renderer, 25, 25, 30, 200);
            SDL_Rect footer = {0, winH-45, 420, 45}; // Expanded width
            SDL_RenderFillRect(renderer, &footer);
            
            SDL_Color c1 = (complexity_level == LEVEL_SKELETAL) ? SDL_Color{0, 255, 100, 255} : SDL_Color{100, 100, 100, 255};
            SDL_Color c2 = (complexity_level == LEVEL_ORGANIC) ? SDL_Color{0, 255, 100, 255} : SDL_Color{100, 100, 100, 255};
            SDL_Color c3 = (complexity_level == LEVEL_LEWIS) ? SDL_Color{0, 255, 100, 255} : SDL_Color{100, 100, 100, 255};
            SDL_Color cR = {255, 100, 100, 255};

            drawText(renderer, font, "[1] Skeletal", 15, by, c1);
            drawText(renderer, font, "[2] Organic", 120, by, c2);
            drawText(renderer, font, "[3] Lewis", 225, by, c3);
            drawText(renderer, font, "[R] Reset", 330, by, cR);
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(16);
    }

    if (font) TTF_CloseFont(font);
    if (fontLarge) TTF_CloseFont(fontLarge);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();

    return 0;
}
