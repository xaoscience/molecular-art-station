#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>

// --- Constants ---
const int SCREEN_WIDTH = 1200;
const int SCREEN_HEIGHT = 800;
const float ATOM_RADIUS = 10.0f;
const float BOND_DISTANCE = 60.0f;
const float MAX_SPEED = 2.0f;

// --- Chemistry Data ---
struct Element {
    std::string symbol;
    int valency; // Max bonds
    SDL_Color color;
};

Element getElement(std::string sym) {
    if (sym == "C") return {"C", 4, {50, 50, 50, 255}};   // Carbon: Dark Grey
    if (sym == "O") return {"O", 2, {200, 50, 50, 255}};  // Oxygen: Red
    if (sym == "N") return {"N", 3, {50, 50, 200, 255}};  // Nitrogen: Blue
    if (sym == "H") return {"H", 1, {200, 200, 200, 255}};// Hydrogen: White
    return {"?", 0, {100, 0, 100, 255}};
}

struct Atom {
    float x, y;
    float vx, vy;
    std::string symbol;
    int max_bonds;
    int current_bonds;
    SDL_Color color;
    int id;
    std::vector<int> bonded_to; // IDs of bonded atoms
};

std::vector<Atom> atoms;
int global_atom_id = 0;

// --- Helper Functions ---

void spawnAtom(float x, float y, std::string sym) {
    Element e = getElement(sym);
    Atom a;
    a.x = x;
    a.y = y;
    // Random velocity
    a.vx = ((float)(rand() % 100) / 50.0f - 1.0f) * MAX_SPEED;
    a.vy = ((float)(rand() % 100) / 50.0f - 1.0f) * MAX_SPEED;
    a.symbol = e.symbol;
    a.max_bonds = e.valency;
    a.current_bonds = 0;
    a.color = e.color;
    a.id = global_atom_id++;
    atoms.push_back(a);
}

float distSq(Atom& a, Atom& b) {
    return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
}

// --- Main Engine ---

int main(int argc, char* args[]) {
    srand(time(NULL));

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    if (TTF_Init() == -1) {
        std::cerr << "TTF could not initialize! TTF_Error: " << TTF_GetError() << std::endl;
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("ChemGen: Molecular Art Station", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    // Load Font (Try to find a system font, or fallback)
    // In a real deployment, we'd ship a font. For now, try common Linux paths.
    TTF_Font* font = TTF_OpenFont("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 16);
    if (!font) {
        // Fallback
        font = TTF_OpenFont("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", 16);
    }
    if (!font) {
        std::cerr << "Warning: Could not load font. Text will not display." << std::endl;
    }

    bool quit = false;
    SDL_Event e;

    // Simulation Loop
    while (!quit) {
        // 1. Input
        int mouseX, mouseY;
        Uint32 mouseState = SDL_GetMouseState(&mouseX, &mouseY);

        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) quit = true;
            if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_ESCAPE) quit = true;
                if (e.key.keysym.sym == SDLK_SPACE) atoms.clear();
            }
        }

        // Mouse Interaction: Spawn atoms
        if (mouseState & SDL_BUTTON(SDL_BUTTON_LEFT)) {
            // Spawn a random organic atom
            int r = rand() % 100;
            std::string sym = "C";
            if (r > 60) sym = "H";
            if (r > 85) sym = "O";
            if (r > 95) sym = "N";
            spawnAtom((float)mouseX + (rand()%20-10), (float)mouseY + (rand()%20-10), sym);
        }

        // 2. Physics & Chemistry
        for (size_t i = 0; i < atoms.size(); i++) {
            // Move
            atoms[i].x += atoms[i].vx;
            atoms[i].y += atoms[i].vy;

            // Bounce off walls
            if (atoms[i].x < 0 || atoms[i].x > SCREEN_WIDTH) atoms[i].vx *= -1;
            if (atoms[i].y < 0 || atoms[i].y > SCREEN_HEIGHT) atoms[i].vy *= -1;

            // Bonding Logic
            // Very naive O(N^2) for now, fine for < 1000 atoms
            if (atoms[i].current_bonds < atoms[i].max_bonds) {
                for (size_t j = i + 1; j < atoms.size(); j++) {
                    if (atoms[j].current_bonds < atoms[j].max_bonds) {
                        // Check if already bonded
                        bool already_bonded = false;
                        for (int id : atoms[i].bonded_to) {
                            if (id == atoms[j].id) already_bonded = true;
                        }
                        
                        if (!already_bonded) {
                            float d2 = distSq(atoms[i], atoms[j]);
                            if (d2 < BOND_DISTANCE * BOND_DISTANCE) {
                                // Form Bond
                                atoms[i].bonded_to.push_back(atoms[j].id);
                                atoms[i].current_bonds++;
                                atoms[j].bonded_to.push_back(atoms[i].id);
                                atoms[j].current_bonds++;
                                
                                // Snap velocities together (inelastic collision)
                                float avg_vx = (atoms[i].vx + atoms[j].vx) / 2.0f;
                                float avg_vy = (atoms[i].vy + atoms[j].vy) / 2.0f;
                                atoms[i].vx = avg_vx;
                                atoms[i].vy = avg_vy;
                                atoms[j].vx = avg_vx;
                                atoms[j].vy = avg_vy;
                            }
                        }
                    }
                }
            }
            
            // Spring forces for bonded atoms (keep them at optimal distance)
            for (int targetId : atoms[i].bonded_to) {
                for (size_t j = 0; j < atoms.size(); j++) {
                    if (atoms[j].id == targetId) {
                        float dx = atoms[j].x - atoms[i].x;
                        float dy = atoms[j].y - atoms[i].y;
                        float dist = sqrt(dx*dx + dy*dy);
                        
                        if (dist > 0) {
                            // Hooke's Law: F = k * (current_dist - optimal_dist)
                            float k = 0.05f; // Spring constant
                            float force = (dist - BOND_DISTANCE) * k;
                            
                            float fx = (dx / dist) * force;
                            float fy = (dy / dist) * force;
                            
                            atoms[i].vx += fx;
                            atoms[i].vy += fy;
                        }
                        break;
                    }
                }
            }

            // Repulsive forces between bonded neighbors (VSEPR approximation)
            // Atoms bonded to the SAME central atom should repel each other
            if (atoms[i].bonded_to.size() > 1) {
                for (size_t b1 = 0; b1 < atoms[i].bonded_to.size(); b1++) {
                    for (size_t b2 = b1 + 1; b2 < atoms[i].bonded_to.size(); b2++) {
                        // Find the two neighbor atoms
                        Atom* n1 = nullptr;
                        Atom* n2 = nullptr;
                        
                        for (size_t k = 0; k < atoms.size(); k++) {
                            if (atoms[k].id == atoms[i].bonded_to[b1]) n1 = &atoms[k];
                            if (atoms[k].id == atoms[i].bonded_to[b2]) n2 = &atoms[k];
                        }
                        
                        if (n1 && n2) {
                            float dx = n1->x - n2->x;
                            float dy = n1->y - n2->y;
                            float distSq = dx*dx + dy*dy;
                            
                            if (distSq > 0 && distSq < (BOND_DISTANCE * 2.5f) * (BOND_DISTANCE * 2.5f)) {
                                float dist = sqrt(distSq);
                                float repulsion = 2.0f / dist; // Inverse distance force
                                
                                float fx = (dx / dist) * repulsion;
                                float fy = (dy / dist) * repulsion;
                                
                                n1->vx += fx;
                                n1->vy += fy;
                                n2->vx -= fx;
                                n2->vy -= fy;
                            }
                        }
                    }
                }
            }

            // Friction / Damping (prevent explosion)
            atoms[i].vx *= 0.98f;
            atoms[i].vy *= 0.98f;
        }

        // 3. Render
        SDL_SetRenderDrawColor(renderer, 20, 20, 25, 255); // Dark background
        SDL_RenderClear(renderer);

        // Draw Bonds
        SDL_SetRenderDrawColor(renderer, 100, 100, 100, 255);
        for (size_t i = 0; i < atoms.size(); i++) {
            for (int targetId : atoms[i].bonded_to) {
                // Find target atom (inefficient search, but simple)
                for (size_t j = 0; j < atoms.size(); j++) {
                    if (atoms[j].id == targetId) {
                        // Draw line
                        SDL_RenderDrawLine(renderer, (int)atoms[i].x, (int)atoms[i].y, (int)atoms[j].x, (int)atoms[j].y);
                        break;
                    }
                }
            }
        }

        // Draw Atoms
        for (auto& a : atoms) {
            // Draw Circle (filled rect for simplicity in low-level SDL without extensions)
            SDL_Rect rect = {(int)(a.x - ATOM_RADIUS), (int)(a.y - ATOM_RADIUS), (int)(ATOM_RADIUS*2), (int)(ATOM_RADIUS*2)};
            SDL_SetRenderDrawColor(renderer, a.color.r, a.color.g, a.color.b, a.color.a);
            SDL_RenderFillRect(renderer, &rect);

            // Draw Symbol
            if (font) {
                SDL_Color textColor = {255, 255, 255};
                if (a.symbol == "H") textColor = {0, 0, 0}; // Black text for Hydrogen

                SDL_Surface* textSurface = TTF_RenderText_Solid(font, a.symbol.c_str(), textColor);
                if (textSurface) {
                    SDL_Texture* textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);
                    SDL_Rect textRect = {(int)(a.x - textSurface->w/2), (int)(a.y - textSurface->h/2), textSurface->w, textSurface->h};
                    SDL_RenderCopy(renderer, textTexture, NULL, &textRect);
                    SDL_FreeSurface(textSurface);
                    SDL_DestroyTexture(textTexture);
                }
            }
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(16); // ~60 FPS
    }

    // Cleanup
    if (font) TTF_CloseFont(font);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();

    return 0;
}
