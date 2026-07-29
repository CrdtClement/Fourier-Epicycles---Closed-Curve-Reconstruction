/*
     ╔═══════════════════════════════════════════════════════════╗
     ║ TIPE 2025-2026, MPI Lycée Louis Thuillier: Cardot Clément ║
     ║  Objective: SDL epicycle animation from Fourier series    ║
     ╚═══════════════════════════════════════════════════════════╝

 ═══════════════════════ Code overview: ════════════════════════════
    STEP 1.
    - Reads a PNG image (passed as argument) using lodepng and extracts
      its contour via the Moore Neighbor Tracing algorithm.
    - Saves the result to sampling.csv as (t, x(t), y(t)) for t in [0;1].

    STEP 2.
    - Computes the 2M+1 complex Fourier coefficients (M given as argument).
    - Outputs fourier_coeffs.csv containing (n, Re(c_n), Im(c_n))
      for all n in [-M ; M].

    STEP 3.
    - Renders an animated epicycle system from the Fourier coefficients
      using an SDL2 graphical window.
    - Displays rendering parameters in real time (circle count, t, scale)
      and allows pausing the animation with the space bar.

 ═════════════════════ Dependencies: ═══════════════════════════════
    - lodepng  : PNG file decoding.
        INCLUDE: lodepng.h and lodepng.c (rename lodepng.cpp if needed).
        Available at https://lodev.org/lodepng/
    - SDL2     : graphical rendering and window management.
        REQUIRES SDL2 (Simple DirectMedia Layer) — https://www.libsdl.org/
        (also SDL2_ttf for on-screen text)
    - math.h   : trigonometric functions.

 ══════════════════════════ Changes from previous version: ══════════
    BEFORE: ./epicycles 0 image.png 50

    The goal was to accept a single image argument and choose M automatically.
    The Riemann sum approximation for Fourier coefficients requires N > M
    (otherwise c_{n+N} = c_n due to periodicity), so M is set to N/10.
    Classic (non-image-based) epicycle mode has been removed.

    NOW:  ./epicycles image.png      => rendering with M := N/10
          ./epicycles image.png M    => rendering with the given M

 ═══════════════════════════ Compilation: ═══════════════════════════
    gcc -o epicycles fourier_epicycles.c lodepng.c -lSDL2 -lSDL2_ttf -lm && ./epicycles image.png
*/


/*
    ╔═══════════════════════════════════════════════════════════╗
    ║                         STEP 1:                           ║
    ║  Read the PNG with lodepng, Moore neighbor tracing algo   ║
    ║           Output format: (t, x(t), y(t))                  ║
    ╚═══════════════════════════════════════════════════════════╝
*/

// ════════════════ include, définitions et types ═══════════════════
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <SDL2/SDL_ttf.h>
#include <math.h>
#include "lodepng.h" 
#include <SDL2/SDL.h>
#define PI 3.14159265358979323846

typedef struct {
    double x, y;
} Point_grille;


// ═════════════════════════ Lecture du PNG ═════════════════════════
unsigned char* decoder_png_rgba(const char* nom_fichier, unsigned* largeur, unsigned* hauteur) {
    unsigned char* image = NULL;
    unsigned erreur = lodepng_decode32_file(&image, largeur, hauteur, nom_fichier);
    if (erreur) {
        printf("Erreur %u: %s\n", erreur, lodepng_error_text(erreur));
        return NULL;
    }
    return image;
}


// ═════════════════════ Traitement de l'image ══════════════════════

unsigned char* convertir_en_N_et_B(unsigned char* rgba, unsigned largeur, unsigned hauteur) {
    unsigned char* gris = malloc(largeur * hauteur);
    for (unsigned y = 0; y < hauteur; y++) {
        for (unsigned x = 0; x < largeur; x++) {
            unsigned indice = 4 * (y * largeur + x);
            unsigned char r = rgba[indice];
            unsigned char g = rgba[indice + 1];
            unsigned char b = rgba[indice + 2];
            unsigned char a = rgba[indice + 3];
            if (a == 0) gris[y*largeur + x] = 255; // Alpha 0 = blanc
            else gris[y*largeur + x] = (r + g + b)/3;
        }
    }
    return gris;
}

void binariser(unsigned char* gris, unsigned largeur, unsigned hauteur, unsigned int seuil) {
    for (unsigned i = 0; i < largeur*hauteur; i++){
        if (gris[i] > seuil)
            gris[i] = 255;
        else
            gris[i] = 0;
    }
}

void inverser_si_necessaire(unsigned char* binaire, unsigned largeur, unsigned hauteur) {
    int noir = 0, blanc = 0;
    for (unsigned i = 0; i < largeur*hauteur; i++) {
        if (binaire[i] == 0)
            noir++;
        else
            blanc++;
    }
    if (blanc < noir) { // la forme est probablement blanche → inverser
        for (unsigned i = 0; i < largeur*hauteur; i++){
            binaire[i] = 255 - binaire[i];
        }
    }
}


// ════════════════ Algorithme du voisinage de Moore ════════════════

int voisins_de_moore(unsigned char* binaire, unsigned largeur, unsigned hauteur, Point_grille** contour, int* max_l, int* max_h) {
    int capacite = 1000, compteur = 0;
    *contour = malloc(capacite * sizeof(Point_grille));

    int x_depart = -1, y_depart = -1;
    for (unsigned y = 1; y < hauteur-1 && x_depart == -1; y++){
        for (unsigned x = 1; x < largeur-1; x++){
            if (binaire[y*largeur + x] == 0){ 
                x_depart = x; 
                y_depart = y; 
                break; 
            }
        }
    }
    if (x_depart == -1) return 0;

    int direction[8][2] = { {0,-1}, {1,-1}, {1,0}, {1,1}, {0,1}, {-1,1}, {-1,0}, {-1,-1} };
    int x_actuel = x_depart, y_actuel = y_depart;
    int x_precedent = x_depart, y_precedent = y_depart - 1;

    int xmin = x_depart, xmax = x_depart, ymin = y_depart, ymax = y_depart;

    do {
        if (compteur >= capacite) {
            capacite *= 2;
            *contour = realloc(*contour, capacite * sizeof(Point_grille));
        }
        (*contour)[compteur++] = (Point_grille){x_actuel, y_actuel};

        if (x_actuel < xmin) xmin = x_actuel;
        if (x_actuel > xmax) xmax = x_actuel;
        if (y_actuel < ymin) ymin = y_actuel;
        if (y_actuel > ymax) ymax = y_actuel;

        int trouve = 0;
        int direction_depart = -1;
        for (int i = 0; i < 8; i++){
            if (x_actuel + direction[i][0] == x_precedent && y_actuel + direction[i][1] == y_precedent){
                direction_depart = (i + 1) % 8; 
                break;
            }
        }
        for (int k = 0; k < 8; k++) {
            int i = (direction_depart + k) % 8;
            int x_suivant = x_actuel + direction[i][0];
            int y_suivant = y_actuel + direction[i][1];
            if (binaire[y_suivant*largeur + x_suivant] == 0) {
                x_precedent = x_actuel;
                y_precedent = y_actuel;
                x_actuel = x_suivant;
                y_actuel = y_suivant;
                trouve = 1;
                break;
            }
        }
        if (!trouve) break;

    } while (x_actuel != x_depart || y_actuel != y_depart);

    *max_l = xmax - xmin;
    *max_h = ymax - ymin;

    return compteur;
}


void calculer_t(Point_grille* contour, int nb_points, double* valeurs_t) {
    for (int i = 0; i < nb_points; i++){
        valeurs_t[i] = (double)i / nb_points;
    }
}


/*  
     ╔═══════════════════════════════════════════════════════════╗
     ║                         ETAPE 2:                          ║
     ║     Calcul coeffs de Fourier par intégrale de Riemann     ║
     ║           Résulat en format (n, Re(c_n), Im(c_n))         ║
     ╚═══════════════════════════════════════════════════════════╝
*/

// ════════════════════ include, def et types ═══════════════════════

typedef struct Complexe {
    double re, im;
} Complexe;

typedef struct Point_parametrique {
    double t, x, y;
} Point_parametrique;

typedef struct Coefficient {
    int n;
    Complexe coeff;
} Coefficient;


// ═════════════════ Primitives sur les complexes ═══════════════════

Complexe complexe_ajout(Complexe a, Complexe b) {
    Complexe r;
    r.re = a.re + b.re;
    r.im = a.im + b.im;
    return r;
}

Complexe complexe_produit(Complexe a, Complexe b) {
    Complexe r;
    r.re = a.re * b.re - a.im * b.im;
    r.im = a.re * b.im + a.im * b.re;
    return r;
}

Complexe complexe_exp(double module, double theta) {
    Complexe r;
    r.re = module * cos(theta);
    r.im = module * sin(theta);
    return r;
}


// ═══════════════════ Calcul des coefficients ══════════════════════

Complexe calculer_un_coeff(Point_parametrique* pts, int N, int n) {
    Complexe somme;
    somme.re = 0.0;
    somme.im = 0.0;
    double T = 1.0; // suppose t dans [0,1]
    double omega = 2*PI/T;

    for (int k = 0; k < N; k++) {
        Complexe ftk;
        ftk.re = pts[k].x;
        ftk.im = pts[k].y;
        double theta = - n * omega * pts[k].t;
        Complexe expterm = complexe_exp(1.0, theta);

        Complexe prod = complexe_produit(ftk, expterm);

        somme.re += prod.re;
        somme.im += prod.im;
    }

    somme.re /= N;
    somme.im /= N;
    return somme;
}

void calculer_tous_les_coeffs(Point_parametrique* pts, int N, int M, Coefficient* coeffs) {
    int idx = 0;
    for (int n = -M; n <= M; n++) {
        Complexe c = calculer_un_coeff(pts, N, n);
        coeffs[idx].n = n;
        coeffs[idx].coeff = c;
        idx++;
    }
}



/*  
     ╔═══════════════════════════════════════════════════════════╗
     ║                         ETAPE 3:                          ║
     ║       Dessin avec SDL, implémentation de l'interface      ║
     ║  def du tableau cercles[] à partir des coeffs de fourier  ║
     ╚═══════════════════════════════════════════════════════════╝
*/


// ════════════════════ include, def et types ═══════════════════════
const int LARGEUR_FENETRE = 800;
const int HAUTEUR_FENETRE = 600;
//const int TAILLE_MAX_CHEMIN = 500; // Voir la ligne disparaitre
const int TAILLE_MAX_CHEMIN = 15000; // Voir la ligne tout le temps
double SPEED = 1;


// ════════════════════════ Struct Cercle ═══════════════════════════
typedef struct {
    double cx, cy;   // centre (double maintenant)
    double rayon;    // rayon (Module du coeff de Fourier)
    double angle;    // phase initiale (argument du coeff de Fourier)
    double vitesse;  // vitesse angulaire (2*PI*n pour le n-ieme coeff, en rad/s)
    int sens;        // sens de rotation : +1 ou -1 (dépend du signe de n)
} Cercle;


// ═════════════════════════ SDL Utils ══════════════════════════════

bool init_SDL(SDL_Window **win, SDL_Renderer **ren) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("Erreur SDL: %s\n", SDL_GetError());
        return false;
    }


    if (TTF_Init() == -1) {
        printf("Erreur TTF_Init: %s\n", TTF_GetError());
        return false;
    }

    *win = SDL_CreateWindow("Epicycles", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                            LARGEUR_FENETRE, HAUTEUR_FENETRE, SDL_WINDOW_SHOWN);
    if (!*win) {
        printf("Erreur creation fenetre: %s\n", SDL_GetError());
        return false;
    }

    *ren = SDL_CreateRenderer(*win, -1, SDL_RENDERER_ACCELERATED);
    if (!*ren) {
        printf("Erreur creation renderer: %s\n", SDL_GetError());
        return false;
    }

    return true;
}


void nettoyage_SDL(SDL_Window *win, SDL_Renderer *ren) {
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    TTF_Quit();
    SDL_Quit();
}

// ══════════════════════════ Dessin ════════════════════════════════

void clear_ecran(SDL_Renderer *ren, int r, int g, int b) {
    SDL_SetRenderDrawColor(ren, r, g, b, 255);
    SDL_RenderClear(ren);
}

void dessiner_cercle(SDL_Renderer *ren, int cx, int cy, int rayon) {
    SDL_SetRenderDrawColor(ren, 90, 90, 90, 255);
    const int NB_SEGMENTS = (2 * PI * rayon) / 5;
    double theta1 = 0;
    double theta2 = (double)1 / NB_SEGMENTS * 2 * PI;
    for (int i = 0; i < NB_SEGMENTS; i++) {
        int x1 = cx + (int)(rayon * cos(theta1));
        int y1 = cy + (int)(rayon * sin(theta1));
        int x2 = cx + (int)(rayon * cos(theta2));
        int y2 = cy + (int)(rayon * sin(theta2));
        SDL_RenderDrawLine(ren, x1, y1, x2, y2);
        theta1 = theta2; 
        theta2 = theta1 + (double)1 / NB_SEGMENTS * 2 * PI;
    }
}

void dessiner_point(SDL_Renderer *ren, Point_grille p, int taille, int r, int g, int b) {
    SDL_SetRenderDrawColor(ren, r, g, b, 100);
    SDL_Rect rect = {(int)(p.x - taille/2), (int)(p.y - taille/2), taille, taille};
    SDL_RenderFillRect(ren, &rect);
}

void dessiner_chemin(SDL_Renderer *ren, Point_grille chemin[], int longueur) {
    SDL_SetRenderDrawColor(ren, 40, 40, 255, 255);
    for (int i = 1; i < longueur; i++) {
        SDL_RenderDrawLine(ren, (int)chemin[i-1].x, (int)chemin[i-1].y, (int)chemin[i].x, (int)chemin[i].y);
    }
}

void afficher_texte(SDL_Renderer *ren, const char *texte, int x, int y, SDL_Color couleur, int taille) {
    TTF_Font *font = TTF_OpenFont("arial.ttf", taille);
    if (!font) {
        printf("Erreur font : %s\n", TTF_GetError());
        return;
    }

    SDL_Surface *surface = TTF_RenderText_Blended(font, texte, couleur);
    SDL_Texture *texture = SDL_CreateTextureFromSurface(ren, surface);

    SDL_Rect rect = {x, y, surface->w, surface->h};
    SDL_RenderCopy(ren, texture, NULL, &rect);

    SDL_FreeSurface(surface);
    SDL_DestroyTexture(texture);
    TTF_CloseFont(font);
}

// ════════════════════════════ Logique ═════════════════════════════

void update_cercle(Cercle *c, double dt) {
    c->angle += c->sens * c->vitesse * SPEED * dt;
    if (c->angle > 2 * PI) c->angle -= 2 * PI;
    if (c->angle < 0) c->angle += 2 * PI;
}

Point_grille point_sur_cercle(Cercle c) {
    Point_grille p;
    p.x = c.cx + c.rayon * cos(c.angle);
    p.y = c.cy + c.rayon * sin(c.angle);
    return p;
}
// ════════════════════ Terminal display ════════════════════════════
const int TEXT_WIDTH = 70;  // inner frame width (adjust if needed)

void print_line(const char *text, int text_width){
    printf(" ║  %-*s ║\n", text_width-3, text); // -3 because "%-*s" counts as 1, offset of 3
}

void print_top_border()
{
    printf("\n\n ╔");
    for (int i = 0; i < TEXT_WIDTH; i++)
        printf("═");
    printf("╗\n");
}

void print_bottom_border()
{
    printf(" ╚");
    for (int i = 0; i < TEXT_WIDTH; i++)
        printf("═");
    printf("╝\n\n\n");
}

void print_empty_line(){
    print_line("", TEXT_WIDTH);
}

void print_title(){
    int space_l = (TEXT_WIDTH - 61)/2;
    int space_r = TEXT_WIDTH - 61 - space_l;

    printf(" ║");
    for(int i = 0; i < space_l; i++){printf(" ");}
    printf("╔");
    for(int i = 0; i < 59; i++){printf("═");}
    printf("╗");
    for(int i = 0; i < space_r; i++){printf(" ");}
    printf("║\n");

    printf(" ║");
    for(int i = 0; i < space_l; i++){printf(" ");}
    printf("║  TIPE 2025-2026, MPI* Lycée Louis Thuillier: Cardot Clément  ║");
    for(int i = 0; i < space_r; i++){printf(" ");}
    printf("║\n");

    printf(" ║");
    for(int i = 0; i < space_l; i++){printf(" ");}
    printf("╚");
    for(int i = 0; i < 59; i++){printf("═");}
    printf("╝");
    for(int i = 0; i < space_r; i++){printf(" ");}
    printf("║\n");
}

void print_usage_error()
{
    print_line("   - image.png : image containing the shape to process", TEXT_WIDTH);
    print_line("                 on a white background", TEXT_WIDTH);
    print_empty_line();
    print_line("    Optional:", TEXT_WIDTH);
    print_line("   - M : Fourier coefficient order (2M+1 coefficients total)", TEXT_WIDTH);
    print_line("         M = N/10 if not specified", TEXT_WIDTH+2);
}


// ══════════════════════ Main loop ═════════════════════════════════

void main_loop(SDL_Renderer *ren, Coefficient *coeffs, int M, float SCALE) {

    int nb_circles = 0;
    Cercle *circles = NULL;

    SPEED = 0.15;
    nb_circles = 2 * M + 1;
    circles = malloc(nb_circles * sizeof(Cercle));
    circles[0] = (Cercle){LARGEUR_FENETRE / 2, HAUTEUR_FENETRE / 2, 0, 0, 0, 1};
    int circle_idx = 1;
    for (int k = 1; k <= M; k++) {
        // -k
        Coefficient coeff_neg = coeffs[M + k];
        circles[circle_idx++] = (Cercle){
            0, 0, SCALE * hypot(coeff_neg.coeff.re, coeff_neg.coeff.im),
            atan2(coeff_neg.coeff.im, coeff_neg.coeff.re), k * 2 * PI, -1
        };
        // +k
        Coefficient coeff_pos = coeffs[M - k];
        circles[circle_idx++] = (Cercle){
            0, 0, SCALE * hypot(coeff_pos.coeff.re, coeff_pos.coeff.im),
            atan2(coeff_pos.coeff.im, coeff_pos.coeff.re), k * 2 * PI, +1
        };
    }

    SDL_Event event;
    bool quit    = false;
    bool paused  = false;
    Uint32 last_time = SDL_GetTicks();
    double t = 0.0;

    Point_grille path[TAILLE_MAX_CHEMIN];
    int path_length = 0;

    // Main drawing loop
    while (!quit) {
        while (SDL_PollEvent(&event)) { // pause / quit events
            if (event.type == SDL_QUIT)
                quit = true;

            if (event.type == SDL_KEYDOWN) {
                switch (event.key.keysym.sym) {
                    case SDLK_ESCAPE:
                        quit = true;
                        break;
                    case SDLK_SPACE:
                        paused = !paused;
                        break;
                }
            }
        }

        // Elapsed time
        Uint32 current_time = SDL_GetTicks();
        double dt = (current_time - last_time) / 1000.0;
        last_time = current_time;

        if (!paused){
            t += dt * SPEED;
            if (t > 1.0){
                t -= 1.0;
            }

            Point_grille current_center = {circles[0].cx, circles[0].cy};

            for (int i = 0; i < nb_circles; i++) {
                circles[i].cx = current_center.x;
                circles[i].cy = current_center.y;
                update_cercle(&circles[i], dt);
                current_center = point_sur_cercle(circles[i]);
            }

            if (path_length < TAILLE_MAX_CHEMIN)
                path[path_length++] = current_center;
            else {
                memmove(path, path + 1, (TAILLE_MAX_CHEMIN - 1) * sizeof(Point_grille));
                path[TAILLE_MAX_CHEMIN - 1] = current_center;
            }
        }

        // Rendering
        clear_ecran(ren, 255, 255, 255);

        Point_grille current_center = {circles[0].cx, circles[0].cy};

        for (int i = 0; i < nb_circles - 1; i++) {
            dessiner_cercle(ren, circles[i].cx, circles[i].cy, circles[i].rayon);
            Point_grille p = point_sur_cercle(circles[i]);
            dessiner_point(ren, p, 1, 180, 0, 0);
            current_center = p;
        }

        dessiner_cercle(ren, circles[nb_circles - 1].cx, circles[nb_circles - 1].cy,
                        circles[nb_circles - 1].rayon);
        Point_grille p = point_sur_cercle(circles[nb_circles - 1]);
        dessiner_point(ren, p, 5, 255, 0, 0);
        current_center = p;

        dessiner_chemin(ren, path, path_length);

        SDL_Color grey = {60, 60, 60, 255};
        SDL_Color green = {0, 255, 100, 255};
        SDL_Color red   = {220, 40, 40, 255};

        char buffer[128];
        snprintf(buffer, sizeof(buffer), "Circles : %d", nb_circles);
        afficher_texte(ren, buffer, 10, 10, grey, 18);

        snprintf(buffer, sizeof(buffer), "SCALE : %.2f", SCALE);
        afficher_texte(ren, buffer, 10, 55, grey, 16);

        snprintf(buffer, sizeof(buffer), "t : %.3f", t);
        afficher_texte(ren, buffer, 10, 75, green, 16);

        if (paused){
            afficher_texte(ren, " PAUSE ", LARGEUR_FENETRE/2 - 45, HAUTEUR_FENETRE - 50, red, 25);
        }

        SDL_RenderPresent(ren);
        SDL_Delay(10);
    }

    free(circles);
}


/*
    ╔═══════════════════════════════════════════════════════════╗
    ║                          MAIN                             ║
    ║                     Shape rendering                       ║
    ╚═══════════════════════════════════════════════════════════╝
*/

int main(int argc, char** argv){
    bool csv = true; // toggle CSV file generation
    char buffer[200]; // used for formatted display lines

    print_top_border();
    print_title();
    print_empty_line();

    if (argc < 2){
        snprintf(buffer, sizeof(buffer), "ERROR: usage: %s image.png", argv[0]);
        print_line(buffer, TEXT_WIDTH);
        print_empty_line();
        print_empty_line();
        print_usage_error();
        print_empty_line();
        print_line("   Exiting.", TEXT_WIDTH);
        print_bottom_border();
        return 1;
    }

    // ══════════════════════════════════════════════════════════
    //          STEP 1: image processing, contour extraction
    // ══════════════════════════════════════════════════════════

    const char* image_file = argv[1];

    snprintf(buffer, sizeof(buffer), "=>  Step 1: Contour extraction from image %s", image_file);
    print_line(buffer, TEXT_WIDTH+2); // +2 for the arrow character width

    unsigned int width, height;
    unsigned char* rgba = decoder_png_rgba(image_file, &width, &height);
    if (!rgba){
        print_line("   ERROR: PNG decoding failed (RGBA is NULL)", TEXT_WIDTH);
        print_bottom_border();
        return 1;
    }

    unsigned char* grey = convertir_en_N_et_B(rgba, width, height);
    free(rgba);

    binariser(grey, width, height, 128);
    inverser_si_necessaire(grey, width, height);

    Point_grille* contour = NULL;

    int max_w, max_h;
    int nb_points = voisins_de_moore(grey, width, height, &contour, &max_w, &max_h);

    double* t_values = malloc(nb_points * sizeof(double));
    calculer_t(contour, nb_points, t_values);

    snprintf(buffer, sizeof(buffer), "    Contour points: %d", nb_points);
    print_line(buffer, TEXT_WIDTH);

    // Write sampling CSV file
    if (csv){
        FILE* sampling_csv = fopen("sampling.csv", "w");
        fprintf(sampling_csv, "t,x,y\n");
        for (int i = 0; i < nb_points; i++){
            fprintf(sampling_csv, "%f,%d,%d\n", t_values[i], (int)(contour[i].x), (int)(contour[i].y));
        }
        fclose(sampling_csv);
        print_line("   sampling.csv written.", TEXT_WIDTH);
    }

    free(grey);

    // ══════════════════════════════════════════════════════════
    //          STEP 2: Fourier coefficient computation
    // ══════════════════════════════════════════════════════════

    int M;
    if (argc > 2){
        M = atoi(argv[2]);
    }
    else{
        M = nb_points / 10; // automatic choice
    }

    if (M > (3 * nb_points) / 4){
        print_line("   WARNING: M too large (> 3/4 of point count, unreliable)", TEXT_WIDTH);
    }

    print_empty_line();
    snprintf(buffer, sizeof(buffer), "=>  Step 2: Fourier coefficient computation (M=%d)", M);
    print_line(buffer, TEXT_WIDTH+2);

    Point_parametrique* points = malloc(nb_points * sizeof(Point_parametrique));
    for (int i = 0; i < nb_points; i++){
        points[i].t = (double)i / nb_points;
        points[i].x = contour[i].x;
        points[i].y = contour[i].y;
    }

    int total = 2*M + 1;
    Coefficient* coeffs = malloc(total * sizeof(Coefficient));
    calculer_tous_les_coeffs(points, nb_points, M, coeffs);

    // Write Fourier coefficients CSV file
    if (csv){
        FILE* fourier_csv = fopen("fourier_coeffs.csv", "w");
        fprintf(fourier_csv, "n,re,im\n");
        for (int i = 0; i < 2*M+1; i++) {
            fprintf(fourier_csv, "%d,%lf,%lf\n", coeffs[i].n, coeffs[i].coeff.re, coeffs[i].coeff.im);
        }
        fclose(fourier_csv);
        print_line("   fourier_coeffs.csv written.", TEXT_WIDTH);
    }

    free(points);
    free(contour);

    // ══════════════════════════════════════════════════════════
    //          STEP 3: Epicycle animation via SDL2
    // ══════════════════════════════════════════════════════════

    print_empty_line();
    print_line("=>  Step 3: Rendering shape with epicycles", TEXT_WIDTH);

    SDL_Window   *win = NULL;
    SDL_Renderer *ren = NULL;
    if (!init_SDL(&win, &ren)){
        print_line("   ERROR: SDL initialisation failed, no window", TEXT_WIDTH);
        print_bottom_border();
        return 1;
    }

    float SCALE = (fmin(LARGEUR_FENETRE, HAUTEUR_FENETRE) * 0.7) / fmax(max_w, max_h);
    snprintf(buffer, sizeof(buffer), "    SCALE = %f", SCALE);
    print_line(buffer, TEXT_WIDTH);
    print_line("    Rendering...", TEXT_WIDTH);
    print_empty_line();

    main_loop(ren, coeffs, M, SCALE);

    nettoyage_SDL(win, ren);

    free(coeffs);
    free(t_values);

    print_line("    Rendering complete.", TEXT_WIDTH);
    print_line("    Exiting.", TEXT_WIDTH);
    print_bottom_border();
    return 0;
}
