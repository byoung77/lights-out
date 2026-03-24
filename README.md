# Lights Out

<p align="center">
  <img src="images/launcher.png" width="45%" />
  <img src="images/gameplay.png" width="45%" />
</p>

An interactive Python implementation of the classic *Lights Out* puzzle, extended with multiple algebraic state spaces and nontrivial topological grids.

---

## Overview

This project implements a fully playable version of *Lights Out* with several enhancements:

- Variable grid sizes from 2×2 up to 15×15
- Multiple state spaces over \( \mathbb{Z}_p \) for \( p \in \{2,3,5,7\} \)
- Support for nontrivial topologies:
  - Simple grid
  - Cylinder
  - Möbius strip
  - Torus
  - Klein bottle
  - Projective plane
- Built-in solver using linear algebra over finite fields
- Graphical user interface with topology visualization

Each game is initialized in a way that guarantees it is solvable.

---

## How It Works

Each board configuration is modeled as a vector over \( \mathbb{Z}_p \), and button presses correspond to adding columns of an adjacency matrix \( A \).

To solve the puzzle, we seek a press vector \( x \) that sends the current board state \( b \) back to the zero state. This means solving the linear system

\[
A x \equiv -b \pmod{p}
\]

where:
- \( b \) is the current board state
- \( x \) is the vector of button presses

A reduced form of the adjacency matrix is computed along with a solver matrix \( A^+ \), allowing solutions to be generated efficiently during gameplay.

If the system has nontrivial nullity, multiple solutions may exist; the program returns one valid solution.

---

## Features

- **Topology-aware neighbor structure**  
  Button interactions depend on the selected surface, including wraparound and twist identifications on spaces such as the torus, Möbius strip, Klein bottle, and projective plane.

- **Guaranteed solvability**  
  Initial states are generated from valid press vectors, ensuring that every puzzle can be completed.

- **Hint system**  
  Generates a valid sequence of button presses to solve the current board.

- **Dynamic UI**  
  - Adjustable board size and state space
  - Visual preview of each topology
  - Responsive grid layout

---

## Requirements

- Python 3.x
- `numpy`
- `Pillow`

Install dependencies with:

```bash
pip install numpy pillow
