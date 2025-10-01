# 3D Lewis Structure and VSEPR Visualizer

This project is a computational chemistry tool written in Python that takes a simple chemical formula, determines its Lewis structure, predicts its 3D molecular geometry using VSEPR (Valence Shell Electron Pair Repulsion) theory, and generates a 3D visualization.

## Features

- **Chemical Formula Parsing**: Parses a chemical formula string (e.g., "H2O", "CH4") into its constituent elements and counts.
- **Central Atom Identification**: Automatically determines the central atom based on electronegativity rules.
- **Lewis Structure Generation**: Calculates total valence electrons and distributes them to form single, double, and triple bonds, as well as lone pairs, to satisfy the octet rule.
- **VSEPR Geometry Prediction**:
  - Determines the steric number for the central atom.
  - Generates an ideal 3D geometry (Linear, Trigonal Planar, Tetrahedral, etc.).
  - Adjusts bond angles to account for the greater repulsion of lone pairs.
- **3D Visualization**: Uses `matplotlib` to create and display an interactive 3D plot of the molecule, representing atoms as spheres and bonds as lines.

## Setup and Usage

Follow these instructions to get the project running on your local machine.

### Prerequisites

- Python 3.8+
- `pip` for package management

### Installation

1.  **Clone the repository:**
    ```bash
    git clone <your-repository-url>
    cd 3D-Lewis
    ```

2.  **Create and activate a virtual environment:**
    ```bash
    # For Windows
    python -m venv venv
    .\venv\Scripts\activate

    # For macOS/Linux
    python3 -m venv venv
    source venv/bin/activate
    ```

3.  **Install dependencies:**
    The project relies on `numpy` for calculations and `matplotlib` for plotting.
    ```bash
    pip install numpy matplotlib
    ```

### Running the Program

To see a demonstration, you can run the `SimpleCompound.py` script directly. The example at the bottom of the file will generate and display the 3D structure for water (H₂O). You can choose the compound by creating a new SimpleCompound object (ex. SimpleCompound("H20")) and calling the displayMolecule().

```bash
python SimpleCompound.py
```

## Future Plans

This project is the foundation for a more advanced molecular modeling tool. Future development will focus on enhancing the accuracy and predictive power of the geometry generation.

1.  **Advanced VSEPR Simulation**: The current model uses a simplified approach to lone pair repulsion. I plan to implement a more physically accurate simulation based on established VSEPR principles where repulsion strengths differ (Lone Pair-Lone Pair > Lone Pair-Bonding Pair > Bonding Pair-Bonding Pair). This will involve a more complex mathematical model to find the minimum energy configuration of the electron domains.

2.  **Machine Learning for Bond Prediction**: While VSEPR provides good estimates, real-world bond angles and lengths are influenced by a multitude of factors. The next major step is to train a machine learning model on a large dataset of known molecular structures (like those from the PubChem database). This model will learn the complex relationships between atomic composition and 3D geometry to predict bond angles and lengths with much higher accuracy than rule-based methods alone.
