# Fourier Epicycles — Closed Curve Reconstruction

> **TIPE 2025-2026 · Cardot Clément · MPI\* · Score: 17.4/20**

## What is a TIPE?

In France, students in *classes préparatoires* (two-year intensive post-baccalaureate programs preparing for the *grandes écoles* competitive exams) must complete a **TIPE** (*Travail d'Initiative Personnelle Encadré* — supervised personal research project). It is evaluated as a standalone oral exam (~15 min presentation + questions) during the national *concours* and accounts for a significant portion of the final score.

The TIPE must be original, scientifically rigorous, and — for the MPI track (Mathematics, Physics, Computer Science) — rooted in at least one of those three disciplines. Each year a broad theme is announced; students must anchor their work to it.
This project was presented as a TIPE oral examination in July 2026 at the IUT de Paris — Université Paris Cité, as part of the MPI* classes préparatoires competitive admissions (concours).
---

## Project overview

**Research question:** To what extent can a closed planar curve, defined by a discrete set of pixels, be faithfully reconstructed by a complex Fourier series — and what are the limits introduced by discretisation, parametrisation, and spectral truncation?

This project implements a full pipeline:

1. **Contour extraction** — Moore neighbor tracing algorithm on a binary PNG image
2. **Fourier analysis** — computation of complex coefficients cₙ via Riemann sum approximation
3. **Epicycle animation** — real-time SDL2 rendering of the rotating circles that trace the curve

The central mathematical result is a proof that the reconstructed curve **converges uniformly** to the original contour (Dirichlet–Jordan theorem), provided the image is sufficiently resolved. The Gibbs phenomenon only arises when the extracted contour is itself discontinuous — a discretisation issue, not a Fourier one.

---

## Results

| Circles | Result |
|---------|--------|
| 9 (sorted by \|cₙ\|) | Rough shape recognisable |
| 101 | Good reconstruction |
| 301 | High-fidelity rendering |

Coefficient selection by modulus significantly improves quality for a fixed number of circles.

---

## Repository structure

```
.
├── src/
│   ├── arial.ttf                 # Arial font file
│   ├── fourier_epicycles.c       # Main C program in English (contour + Fourier + SDL2)
│   ├── fourier_epicycles_fr.c    # Main C program in French (contour + Fourier + SDL2)
│   ├── lodepng.h                 # LodePNG header (PNG decoder)
│   ├── lodepng.c                 # LodePNG implementation
│   └── plot_contour.py           # Python script for contour visualisation
├── docs/
│   ├── presentation.pdf          # Oral presentation slides
│   └── MCOT.pdf                  # Research objectives document (MCOT in French)
├── examples/
│   └── *.png                     # Sample input images
└── README.md
```

---

## Dependencies

| Dependency | Purpose |
|------------|---------|
| [SDL2](https://www.libsdl.org/) | Real-time graphics rendering |
| [SDL2_ttf](https://github.com/libsdl-org/SDL_ttf) | On-screen text display |
| [LodePNG](https://lodev.org/lodepng/) | PNG decoding (included in `src/`) |
| `libm` | Math functions |
| Python + matplotlib + numpy | Contour plotting (optional) |

### Installing SDL2 (Ubuntu/Debian)

```bash
sudo apt install libsdl2-dev libsdl2-ttf-dev
```

### Installing SDL2 (macOS)

```bash
brew install sdl2 sdl2_ttf
```

---

## Build & run

```bash
gcc -o epicycles src/fourier_epicycles.c src/lodepng.c -lSDL2 -lSDL2_ttf -lm
```

```bash
# Automatic mode: M = N/10 (N = number of contour points)
./epicycles examples/pi.png

# Manual mode: specify the number of Fourier coefficients (2M+1 total)
./epicycles examples/pi.png 50
```

**Controls:**
- `Space` — pause/resume the animation
- `Escape` — quit

---

## How it works

### 1 — Contour extraction (Moore neighbor tracing)

The input image is converted to greyscale, binarised, and then traced using the **Moore neighbor algorithm**: starting from the first black pixel found, the algorithm walks along the boundary by checking the 8-connected neighborhood at each step, always entering from the direction of the previous pixel. This guarantees a complete, ordered, closed contour.

The result is a sequence of N points P_k = (x_k, y_k) with t_k = k/N ∈ [0, 1].

### 2 — Fourier coefficient computation

Each point is encoded as a complex number: γ(t_k) = x_k + i·y_k.

The interpolated curve γ is piecewise affine between consecutive points, so γ ∈ C⁰([0,1]) ∩ C¹_piecewise([0,1]).

The Fourier coefficients are computed by Riemann sum:

```
cₙ = (1/N) · Σ γ(t_k) · exp(−2πi·n·t_k)
```

**Sampling constraint:** the discrete sum satisfies cₙ = cₙ₊ₙ (periodicity), so M must satisfy M ≤ N to avoid aliasing.

**Convergence:** by the **Dirichlet–Jordan theorem**, the partial sums S_M(γ)(t) converge uniformly to γ as M → ∞.

### 3 — Epicycle rendering (SDL2)

Each Fourier coefficient cₙ = |cₙ|·e^(iθₙ) defines a circle:
- **radius** = |cₙ| (geometric importance of the frequency)
- **initial angle** = arg(cₙ)
- **angular velocity** = 2π·|n| rad/s, rotating clockwise for n < 0

The circles are chained: the centre of each circle sits on the tip of the previous one. The tip of the last circle traces the reconstructed curve.

---

## Identified limits

- **Gibbs phenomenon** — appears only when the extracted contour is discontinuous (image under-resolved, typically < ~100 px). Does not occur for smooth or piecewise-smooth curves.
- **Uniform parametrisation** — t_k = k/N may introduce artefacts if the contour has highly non-uniform point density. Arc-length parametrisation would be more accurate.
- **Single contour** — the algorithm extracts one boundary per run; images with holes or multiple shapes are not handled.

---

## Academic context

- **Programme:** MPI\* (Mathematics, Physics, Computer Science), 2nd year (*5/2*)
- **Theme:** *Cycles et Boucles* (Cycles and Loops) — 2025-2026
- **Thematic positioning:** Computer Science (practical) · Mathematics (Applied)
- **Keywords:** Fourier series · Moore neighborhood · Discrete contour · Planar parametric curve · Epicycles

---

## References

1. Seo et al. — *Fast Contour-Tracing Algorithm Based on a Pixel-Following Method* (2016)
2. David Doria — *Moore-Neighbor Tracing*, Insight Journal
3. Stein & Shakarchi — *Fourier Analysis: An Introduction*, Princeton University Press
4. Hanson — *The Mathematical Power of Epicyclical Astronomy*, Isis (1960)
5. Mayer, Khairy, Howard — *Drawing an elephant with four complex parameters*, AJP (2010)
6. [SDL2 Wiki](https://wiki.libsdl.org/SDL2/FrontPage)
