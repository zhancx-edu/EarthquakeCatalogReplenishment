A-BETA: Fully Automated Earthquake Catalog Replenishment Tool
Overview
This repository contains the Python implementation and data for A-BETA (Auto-Biscale Empirical Transformation Application), a fully automated tool for replenishing incomplete earthquake catalogs, as described in the paper submitted to Seismological Research Letters (SRL, under review). The tool enhances seismological analysis by addressing missing early aftershock data through grid-based density analysis and flexible missing region identification.
Paper Details

Title: A-BETA: Fully Automated Earthquake Catalog Replenishment Tool
Authors: Chengxiang Zhan, Jiancang Zhuang, Stephen Wu
Abstract:
Incomplete earthquake catalogs, especially for early aftershocks, pose challenges for accurate seismological analysis. This study introduces Auto-Biscale Empirical Transformation Application (A-BETA), a fully automated earthquake data replenishment method that builds upon previous semi-automated approaches by integrating grid-based density analysis with flexible missing region identification. The method supports user-defined polygons and hybrid strategies (union or intersection) to incorporate expert knowledge. We applied A-BETA to the Wenchuan ($M_w$ 7.9, 2008) and Noto ($M_w$ 7.5, 2024) earthquake sequences, using magnitude thresholds of $M \geq 3.0$ and $M \geq 2.0$, respectively. After 100 automated iterations, the method replenished 43.0% and 39.8% of the final catalogs. These results substantially improve catalog completeness, enabling more robust aftershock modeling and statistical analysis. A-BETA is implemented in Python and is available as open-source software on GitHub.



Repository Contents

src/: Source code for the A-BETA tool, written in Python.
data/: Sample datasets for the Wenchuan (2008, $M_w$ 7.9) and Noto (2024, $M_w$ 7.5) earthquake sequences.
examples/: Example scripts demonstrating how to use A-BETA with provided datasets.
docs/: Additional documentation (if applicable, e.g., user guide or API reference).
.gitignore: Configuration to exclude unnecessary files (e.g., __pycache__, virtual environments).
requirements.txt: Python dependencies required to run A-BETA.
LICENSE: License file (e.g., MIT, to be specified).

Installation

Clone the repository:git clone https://github.com/<your-username>/<your-repo-name>.git
cd <your-repo-name>


Create and activate a virtual environment (optional but recommended):python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate


Install dependencies:pip install -r requirements.txt


Verify installation by running an example script:python examples/wenchuan_example.py



Usage
A-BETA supports user-defined polygons and hybrid strategies (union or intersection) for missing region identification. To apply A-BETA to your own earthquake catalog:

Prepare your dataset in the format specified in docs/data_format.md (e.g., CSV with columns for time, latitude, longitude, magnitude).
Configure parameters (e.g., magnitude threshold, iteration count) in the script or via command-line arguments.
Run the tool:python src/abeta.py --input data/your_catalog.csv --output results/replenished_catalog.csv



See examples/ for detailed usage with the Wenchuan and Noto datasets.
Results
When applied to:

Wenchuan (2008): Replenished 43.0% of the catalog ($M \geq 3.0$) after 100 iterations.
Noto (2024): Replenished 39.8% of the catalog ($M \geq 2.0$) after 100 iterations.

These results enhance catalog completeness for aftershock modeling and statistical analysis.
Contributing
Contributions are welcome! Please:

Fork the repository.
Create a feature branch (git checkout -b feature/your-feature).
Commit changes (git commit -m "Add your feature").
Push to the branch (git push origin feature/your-feature).
Open a Pull Request.

License
This project is licensed under the MIT License (or specify another license as appropriate).
Citation
If you use A-BETA in your research, please cite our paper (pending publication in Seismological Research Letters):

Zhan, C., Zhuang, J., & Wu, S. (2025). A-BETA: Fully Automated Earthquake Catalog Replenishment Tool. Seismological Research Letters (under review).

Contact
For questions or issues, please open an issue on GitHub or contact the authors:

Chengxiang Zhan: [email TBD]
Jiancang Zhuang: [email TBD]
Stephen Wu: [email TBD]

