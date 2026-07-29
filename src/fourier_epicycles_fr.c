/*  
     ╔═══════════════════════════════════════════════════════════╗
     ║ TIPE 2025-2026,MPI* Lycée Louis Thuillier: Cardot Clément ║
     ║  Objectif: Animation SDL d'épicycles à partir de Fourier  ║
     ╚═══════════════════════════════════════════════════════════╝

 ═══════════════════════ Objectif de ce code: ═══════════════════════
    ETAPE 1.
    -Permet de lire une image PNG (passé en paramètre) grâce à lodepng et en extraire le contour par l'algorithme Moore Neighbor Tracing.
    -Enregistre le résultat dans le fichier echantillonage.csv (t, x(t), y(t)) pour t dans [0;1].
        
    ETAPE 2.
    -calcule les 2M + 1 coefficients de Fourier complexes associés (M donné en paramètre). 
    -En sortie, le programme génère le fichier coeffs_fourier.csv contenant (n, Re(c_n), Im(c_n)) pour tout n dans [-M ; M].

    ETAPE 3.
    -Permet de tracer des épicycles animés à partir de coefficients de Fourier, avec une interface graphique réalisée avec SDL2.
    -Affiche en temps réel les paramètres du tracé (nombre de cercles, t, scale) et permet de mettre en pause l'animation avec la barre d'espace.

 ═════════════════════ Bibliothèques utilisées: ═════════════════════
    - lodepng pour la lecture des fichiers PNG.
        JOINDRE: lodepng.h et lodepng.c (lodepng.cpp à renommer), fichiers sur https://lodev.org/lodepng/
    - SDL2 pour le rendu graphique et la gestion de la fenêtre.
        NECESSITE SDL2 (Simple DirectMedia Layer) - https://www.libsdl.org/
        (aussi SDL2_ttf pour le texte)
    - math.h pour les fonctions trigonométriques.

 ══════════════════════════ Modifications: ══════════════════════════
    AVANT: ./epicycles 0 image.png 50 
    
    L'objectif est donc de prendre seulement une image en argument, contenant le contour à tracer et de faire le choix de M automatiquement.
    L'approximation des intégrales pour le calcul des coefficients de Fourier oblige à respecter N > M (sinon c_{n+N} = c_n), donc je fixe M:=N/10.
    Plus possible de tracer les épicycles classiques.

    MAINTENANT: ./epicycles image.png     => tracé avec M:=N/10
                ./epicycles image.png M   => tracé avec M donné en paramètres

 ═══════════════════════════ Compilation: ═══════════════════════════
    gcc -o epicycles fourier_epicycles_fr.c lodepng.c -lSDL2 -lSDL2_ttf -lm && ./Dessin image.png
*/


/*  
     ╔═══════════════════════════════════════════════════════════╗
     ║                         ETAPE 1:                          ║
     ║ Lire le .png avec lodepng, algorithe des voisins de Moore ║
     ║              Résulat en format (t, x(t), y(t))            ║
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



// ════════════════════ Affichage du terminal ═══════════════════════
const int LARGEUR_TEXTE = 70;  // largeur intérieure du cadre (à ajuster si besoin)

void afficher_ligne(const char *texte, int largeur_texte){
    printf(" ║  %-*s ║\n", largeur_texte-3, texte); //-3 car "%-*s" est considéré comme 1 donc decalage de 3...
}

void afficher_bord_haut()
{
    printf("\n\n ╔");
    for (int i = 0; i < LARGEUR_TEXTE; i++)
        printf("═");
    printf("╗\n");
}

void afficher_bord_bas()
{
    printf(" ╚");
    for (int i = 0; i < LARGEUR_TEXTE; i++)
        printf("═");
    printf("╝\n\n\n");
}
void afficher_ligne_vide(){
    afficher_ligne("", LARGEUR_TEXTE);
}
void afficher_titre(){
    int espace_g = (LARGEUR_TEXTE - 61)/2;
    int espace_d = LARGEUR_TEXTE - 61 -espace_g;
    printf(" ║");
    for(int i = 0; i < espace_g; i++){printf(" ");}
    printf("╔");
    for(int i = 0; i < 59; i++){printf("═");}
    printf("╗");
    for(int i = 0; i < espace_d; i++){printf(" ");}
    printf("║\n");

    printf(" ║");
    for(int i = 0; i < espace_g; i++){printf(" ");}
    printf("║ TIPE 2025-2026,MPI* Lycée Louis Thuillier:Cardot Clément ║");
    for(int i = 0; i < espace_d; i++){printf(" ");}
    printf("║\n");

    printf(" ║");
    for(int i = 0; i < espace_g; i++){printf(" ");}
    printf("╚");
    for(int i = 0; i < 59; i++){printf("═");}
    printf("╝");
    for(int i = 0; i < espace_d; i++){printf(" ");}
    printf("║\n");
}

void afficher_erreur_usage()
{
    afficher_ligne("   - image.png : image contenant la forme a traiter", LARGEUR_TEXTE);
    afficher_ligne("                 sur fond blanc", LARGEUR_TEXTE);
    afficher_ligne_vide();
    afficher_ligne("    Facultatif: ", LARGEUR_TEXTE);
    afficher_ligne("   - M : ordre des coefficients de Fourier (2M+1 coeffs)", LARGEUR_TEXTE);
    afficher_ligne("         M = N/10 si pas spécifié", LARGEUR_TEXTE+2);
}


// ══════════════════════ Boucle principale ═════════════════════════

void main_loop(SDL_Renderer *ren, Coefficient *coeffs, int M, float SCALE) {
    
    int nb_cercles = 0;
    Cercle *cercles = NULL;

    SPEED = 0.15;
    nb_cercles = 2 * M + 1;
    cercles = malloc(nb_cercles * sizeof(Cercle));
    cercles[0] = (Cercle){LARGEUR_FENETRE / 2, HAUTEUR_FENETRE / 2, 0, 0, 0, 1};
    int num_cercle = 1;
    for (int k = 1; k <= M; k++) {
        // -k
        Coefficient coeff_neg = coeffs[M + k];
        cercles[num_cercle++] = (Cercle){
            0, 0, SCALE * hypot(coeff_neg.coeff.re, coeff_neg.coeff.im),atan2(coeff_neg.coeff.im, coeff_neg.coeff.re),k * 2 * PI, -1
        };
        // +k
        Coefficient coeff_pos = coeffs[M - k];
        cercles[num_cercle++] = (Cercle){
            0, 0, SCALE * hypot(coeff_pos.coeff.re, coeff_pos.coeff.im),atan2(coeff_pos.coeff.im, coeff_pos.coeff.re),k * 2 * PI, +1
        };
    }

    SDL_Event epicycles;
    bool quitter = false;
    bool en_pause = false;
    Uint32 last_time = SDL_GetTicks();
    double t = 0.0;

    Point_grille chemin[TAILLE_MAX_CHEMIN];
    int longueur_chemin = 0;

    //  Boucle principal de traçage 
    while (!quitter) {
        while (SDL_PollEvent(&epicycles)) { //Action pause
            if (epicycles.type == SDL_QUIT)
                quitter = true;

            if (epicycles.type == SDL_KEYDOWN) {
                switch (epicycles.key.keysym.sym) {
                    case SDLK_ESCAPE:
                        quitter = true;
                        break;
                    case SDLK_SPACE:
                        en_pause = !en_pause; 
                        break;
                }
            }
        }

        //Calcul temps ecoulé
        Uint32 current_time = SDL_GetTicks();
        double dt = (current_time - last_time) / 1000.0;
        last_time = current_time;

        if(!en_pause){
            t += dt * SPEED;
            if(t > 1.0){
                t -= 1.0;
            }

            Point_grille centre_courant = {cercles[0].cx, cercles[0].cy};

            for (int i = 0; i < nb_cercles; i++) {
                cercles[i].cx = centre_courant.x;
                cercles[i].cy = centre_courant.y;
                update_cercle(&cercles[i], dt);
                centre_courant = point_sur_cercle(cercles[i]);
            }

            if (longueur_chemin < TAILLE_MAX_CHEMIN)
                chemin[longueur_chemin++] = centre_courant;
            else {
                memmove(chemin, chemin + 1, (TAILLE_MAX_CHEMIN - 1) * sizeof(Point_grille));
                chemin[TAILLE_MAX_CHEMIN - 1] = centre_courant;
            }
        }

        // Rendu graphique
        clear_ecran(ren, 255, 255, 255);

        Point_grille centre_courant = {cercles[0].cx, cercles[0].cy};

        for (int i = 0; i < nb_cercles - 1; i++) {
            dessiner_cercle(ren, cercles[i].cx, cercles[i].cy, cercles[i].rayon);
            Point_grille p = point_sur_cercle(cercles[i]);
            dessiner_point(ren, p, 4, 255, 0, 128);
            centre_courant = p;
        }

        dessiner_cercle(ren, cercles[nb_cercles - 1].cx, cercles[nb_cercles - 1].cy,
                        cercles[nb_cercles - 1].rayon);
        Point_grille p = point_sur_cercle(cercles[nb_cercles - 1]);
        dessiner_point(ren, p, 7, 255, 0, 128);
        centre_courant = p;

        dessiner_chemin(ren, chemin, longueur_chemin);

        SDL_Color gris = {60, 60, 60, 255};
        SDL_Color vert = {0, 255, 100, 255};
        SDL_Color rouge = {220, 40, 40, 255};

        char buffer[128];
        snprintf(buffer, sizeof(buffer), "Cercles : %d", nb_cercles);
        afficher_texte(ren, buffer, 10, 10, gris, 16);

        snprintf(buffer, sizeof(buffer), "SCALE : %.2f", SCALE);
        afficher_texte(ren, buffer, 10, 55, gris, 16);

        snprintf(buffer, sizeof(buffer), "t : %.3f", t);
        afficher_texte(ren, buffer, 10, 75, vert, 16);

        if(en_pause){
            afficher_texte(ren, " PAUSE ", LARGEUR_FENETRE/2 - 45, HAUTEUR_FENETRE - 50, rouge, 25);
        }

        SDL_RenderPresent(ren);
        SDL_Delay(10);
    }

    free(cercles);
}




/*  
     ╔═══════════════════════════════════════════════════════════╗
     ║                           MAIN                            ║
     ║                    Traçage de la forme                    ║
     ║                                                           ║
     ╚═══════════════════════════════════════════════════════════╝

*/

int main(int argc, char** argv){
    bool csv = false; // Permet de choisir si on veut générer les fichiers csv ou non
    char buffer[200]; //Gerer l'affichage des espaces

    afficher_bord_haut();
    afficher_titre();
    afficher_ligne_vide();

    if (argc < 2){
        snprintf(buffer, sizeof(buffer), "ERREUR: usage: %s image.png", argv[0]);
        afficher_ligne(buffer, LARGEUR_TEXTE);
        afficher_ligne_vide();
        afficher_ligne_vide();
        afficher_erreur_usage();
        afficher_ligne_vide();
        afficher_ligne("   Fermeture du programme.", LARGEUR_TEXTE);
        afficher_bord_bas();
        return 1;
    }

    // ══════════════════════════════════════════════════════════════════
    //         ETAPE 1: on traite l'image, calcul du contour
    // ══════════════════════════════════════════════════════════════════

    const char* image_file = argv[1];

    snprintf(buffer, sizeof(buffer), "⇾  Étape 1 : Extraction du contour de l'image %s", image_file);
    afficher_ligne(buffer, LARGEUR_TEXTE+3); //+3 car la flèche est un caractère qui compte pour 2, %s aussi et l'accent aussi...
    
    unsigned int largeur, hauteur;
    unsigned char* rgba = decoder_png_rgba(image_file, &largeur, &hauteur);
    if (!rgba){
        afficher_ligne("   ERREUR PNG, RGBA FALSE", LARGEUR_TEXTE);
        afficher_bord_bas();
        return 1;
    } 

    unsigned char* gris = convertir_en_N_et_B(rgba, largeur, hauteur);
    free(rgba);

    binariser(gris, largeur, hauteur, 128);
    inverser_si_necessaire(gris, largeur, hauteur);

    Point_grille* contour = NULL;

    int max_l, max_h;
    int nb_points = voisins_de_moore(gris, largeur, hauteur, &contour, &max_l, &max_h);

    double* valeurs_t = malloc(nb_points * sizeof(double));
    calculer_t(contour, nb_points, valeurs_t);
    snprintf(buffer, sizeof(buffer), "    Nombre de points du contour : %d", nb_points);
    afficher_ligne(buffer, LARGEUR_TEXTE);

    // Écrire le fichier echantillonage.csv
    if(csv){
        FILE* echantillonage_csv = fopen("echantillonage.csv", "w");
        fprintf(echantillonage_csv, "t,x,y\n");
        for (int i = 0; i < nb_points; i++){
            fprintf(echantillonage_csv, "%f,%d,%d\n", valeurs_t[i], (int)(contour[i].x), (int)(contour[i].y));
        }
        fclose(echantillonage_csv);
        afficher_ligne("   Fichier echantillonage.csv écrit.", LARGEUR_TEXTE+1); // "é" compte pour 2 (Pourquoi ?! aucune idée...)
    }

    free(gris);

    // ══════════════════════════════════════════════════════════════════
    //            ETAPE 2: Calcul des coefficients de Fourier
    // ══════════════════════════════════════════════════════════════════
    
    int M;
    if (argc > 2){
        M = atoi(argv[2]);
    }
    else{
        M = nb_points / 10;   // choix automatique
    }
    afficher_ligne_vide();
    snprintf(buffer, sizeof(buffer), "⇾  Étape 2 : Calcul des coefficients de Fourier (M=%d)", M);
    afficher_ligne(buffer, LARGEUR_TEXTE+3); //idem, fleche et accent


    Point_parametrique* points = malloc(nb_points * sizeof(Point_parametrique));
    for (int i = 0; i < nb_points; i++){
        points[i].t = (double)i / nb_points;
        points[i].x = contour[i].x;
        points[i].y = contour[i].y;
    }
    int total = 2*M + 1;
    Coefficient* coeffs = malloc(total * sizeof(Coefficient));
    calculer_tous_les_coeffs(points, nb_points, M, coeffs);

    // Écrire le fichier coeffs_fourier.csv
    if(csv){
        FILE* coeffs_fourier_csv = fopen("coeffs_fourier.csv", "w");
        fprintf(coeffs_fourier_csv, "n,re,im\n");
        for (int i = 0; i < 2*M+1; i++) {
            fprintf(coeffs_fourier_csv, "%d,%lf,%lf\n", coeffs[i].n, coeffs[i].coeff.re, coeffs[i].coeff.im);
        }
        fclose(coeffs_fourier_csv);
        afficher_ligne("    Fichier coeffs_fourier.csv écrit.", LARGEUR_TEXTE+1); // "é" compte pour 2
    }

    free(points);
    free(contour);

    // ══════════════════════════════════════════════════════════════════
    //        ETAPE 3: Dessin des épicycles à l'aide des coeffs
    // ══════════════════════════════════════════════════════════════════

    afficher_ligne_vide();
    afficher_ligne("⇾  Étape 3 : Tracé de la forme avec des epicycles", LARGEUR_TEXTE+4); //fleches et accents...

    SDL_Window *win = NULL;
    SDL_Renderer *ren = NULL;
    if (!init_SDL(&win, &ren)){
        afficher_ligne("   ERREUR SDL, PAS DE FENETRE", LARGEUR_TEXTE);
        afficher_bord_bas();
        return 1;
    }

    float SCALE = (fmin(LARGEUR_FENETRE, HAUTEUR_FENETRE) * 0.7) / fmax(max_l, max_h);
    snprintf(buffer, sizeof(buffer), "    SCALE = %f", SCALE);
    afficher_ligne(buffer, LARGEUR_TEXTE);
    afficher_ligne("    Affichage en cours...", LARGEUR_TEXTE);
    afficher_ligne_vide();
    main_loop(ren, coeffs, M, SCALE);

    nettoyage_SDL(win, ren);
    
    free(coeffs);
    free(valeurs_t);

    afficher_ligne("    Fin affichage.", LARGEUR_TEXTE);
    afficher_ligne("    Fermeture du programme.", LARGEUR_TEXTE);
    afficher_bord_bas();
    return 0;
}
