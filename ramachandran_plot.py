import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from fpdf import FPDF, XPos, YPos

# Función para calcular el ángulo de torsión entre cuatro átomos
def torsion_angle(a, b, c, d):
    # Calcula vectores entre los puntos
    u1 = np.array(b) - np.array(a)
    u2 = np.array(c) - np.array(b)
    u3 = np.array(d) - np.array(c)

    # Calcula productos cruzados y puntos necesarios para el ángulo
    e = (np.linalg.norm(u2) * u1)
    f = np.cross(u2, u3)
    g = np.cross(u1, u2)

    # Calcula componentes del ángulo
    y = np.dot(e, f)
    x = np.dot(g, f)

    # Calcula y retorna el ángulo en grados
    ang = np.arctan2(y, x)
    return np.degrees(ang)

# Función para obtener las coordenadas de un átomo específico
def get_atom_coords(df, chain, residue, atom):
    try:
        atom_row = df[(df["chain"] == chain) & (df["residue num"] == residue) & (df["atom"] == atom)]
        return atom_row[['x', 'y', 'z']].values[0]
    except IndexError:
        return None

# Función para clasificar residuos de aminoácidos
def classify_residue(residue):
    if residue in polar:
        return "Polar"
    elif residue in polar_pos:
        return "Polar positive charged"
    elif residue in polar_neg:
        return "Polar negative charged"
    elif residue in nonpolar_ali:
        return "Non polar aliphatic"
    elif residue in nonpolar_aro:
        return "Non polar aromatic"
    else:
        pass

# Nombre del archivo PDB y otras configuraciones
pdb = "1ewq.pdb"
prot_pdb = pdb.split(".")[0].upper()
data = []
angles = []

# Listas de aminoácidos y sus clasificaciones
aminoacids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
              "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
polar = ["ASN", "CYS", "GLN", "SER", "THR"]
polar_pos = ["ARG", "HIS", "LYS"]
polar_neg = ["ASP", "GLU"]
nonpolar_ali = ["ALA", "ILE", "GLY", "LEU", "MET", "PRO", "VAL"]
nonpolar_aro = ["PHE", "TYR", "TRP"]

# Leer el archivo PDB y extraer datos relevantes
with open(pdb) as file:
    for line in file:
        if line.startswith("ATOM"):
            atom = line[12:16].strip()
            residue = line[17:20].strip()
            chain = line[21].strip()
            res_num = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            if residue in aminoacids:
                data.append([atom, residue, chain, res_num, x, y, z])

# Crear un DataFrame a partir de los datos extraídos
titles = ["atom", "residue", "chain", "residue num", "x", "y", "z"]
df = pd.DataFrame(data, columns=titles)

# Calcular ángulos phi y psi para cada residuo en cada cadena
for chain in df["chain"].unique():
    for residue in df[df["chain"] == chain]["residue num"].unique():
        residue_name = df[(df["chain"] == chain) & (df["residue num"] == residue)]["residue"].values[0]
        C_ant = get_atom_coords(df, chain, residue - 1, "C")
        N = get_atom_coords(df, chain, residue, "N")
        CA = get_atom_coords(df, chain, residue, "CA")
        C = get_atom_coords(df, chain, residue, "C")
        N_next = get_atom_coords(df, chain, residue + 1, "N")

        phi = float()
        psi = float()

        if C_ant is not None and N is not None and CA is not None and C is not None:
            phi = round(torsion_angle(C_ant, N, CA, C), 2)

        if N is not None and CA is not None and C is not None and N_next is not None:
            psi = round(torsion_angle(N, CA, C, N_next), 2)

        if phi != 0.00 and psi != 0.00:
            angles.append((chain, residue, residue_name, phi, psi))

# Crear un DataFrame a partir de los ángulos calculados
angles_df = pd.DataFrame(angles, columns=["chain", "residue num", "residue", "phi", "psi"])

# Clasificar los residuos
angles_df["classification"] = angles_df["residue"].apply(classify_residue)

# Filtrar residuos específicos (prolina y glicina)
proline_angles_df = angles_df.loc[angles_df['residue'] == 'PRO']
glycine_angles_df = angles_df.loc[angles_df['residue'] == 'GLY']

# Generar y guardar el gráfico de Ramachandran para glicina
glycine_plot_name = f"{prot_pdb}_ramachandran_plot_glycine.png"
plt.figure(figsize=(8, 6))
sns.scatterplot(x="phi", y="psi", data=glycine_angles_df, s=15, edgecolor="green", color="lightgreen")
sns.kdeplot(x="phi", y="psi", data=glycine_angles_df, levels=3, color="darkgreen", fill=True, alpha=0.3, bw_adjust=0.25)
plt.xlabel("Phi (φ)")
plt.ylabel("Psi (ψ)")
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.title("Ramachandran plot for glycine residues")
plt.grid(True)
plt.savefig(glycine_plot_name)

# Generar y guardar el gráfico de Ramachandran para prolina
proline_plot_name = f"{prot_pdb}_ramachandran_plot_proline.png"
plt.figure(figsize=(8, 6))
sns.scatterplot(x="phi", y="psi", data=proline_angles_df, s=15, edgecolor="red", color="tomato")
sns.kdeplot(x="phi", y="psi", data=proline_angles_df, levels=3, color="darkred", fill=True, alpha=0.3, bw_adjust=0.25)
plt.xlabel("Phi (φ)")
plt.ylabel("Psi (ψ)")
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.title("Ramachandran plot for proline residues")
plt.grid(True)
plt.savefig(proline_plot_name)

# Generar y guardar el gráfico de Ramachandran general
plot_name = f"{prot_pdb}_ramachandran_plot.png"
plt.figure(figsize=(8, 6))
sns.scatterplot(x="phi", y="psi", data=angles_df, s=15, edgecolor="blue", color="lightblue")
sns.kdeplot(x="phi", y="psi", data=angles_df, levels=4, color="darkblue", fill=True, alpha=0.3)
plt.xlabel("Phi (φ)")
plt.ylabel("Psi (ψ)")
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.title("Ramachandran plot")
plt.grid(True)
plt.savefig(plot_name)

# Generar y guardar el histograma de frecuencias de residuos
hist_plot_name = f"{prot_pdb}_residue_histogram.png"
plt.figure(figsize=(8, 6))
ax = sns.histplot(angles_df["residue"], color="lightblue")
for p in ax.patches:
    ax.text(p.get_x() + p.get_width() / 2., -0.05,
            f'{int(p.get_height())}', ha='center', va='center', fontsize=8,
            bbox=dict(facecolor='white', alpha=0.9))
plt.xlabel("Residue")
plt.ylabel("Repetitions")
plt.title("Frecuency of residues in the Ramachandran plot")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(hist_plot_name)

# Generar y guardar el gráfico de clasificación de residuos
classification_name = f"{prot_pdb}_classification.png"
plt.figure(figsize=(8, 6))
sns.countplot(x="classification", data=angles_df, palette="Paired", order=angles_df["classification"].value_counts().index)
plt.xlabel("Classification")
plt.ylabel("Repetitions")
plt.title("Classification of the residues in the Ramachandran plot")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(classification_name)

# Clase para crear el PDF con los gráficos
class PDF(FPDF):
    def header(self):
        self.set_font('Helvetica', 'B', 12)
        self.cell(0, 10, prot_pdb, new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')

    def footer(self):
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}', new_x=XPos.RIGHT, new_y=YPos.TOP, align='C')

# Crear el PDF y añadir los gráficos
pdf = PDF()
pdf.add_page()
pdf.image(plot_name, x=5, y=19, w=205)
pdf.image(proline_plot_name, x=10, y=180, w=100)
pdf.image(glycine_plot_name, x=105, y=180, w=100)
pdf.add_page()
pdf.image(hist_plot_name, x=10, y=19, w=180)
pdf.image(classification_name, x=10, y=160, w=180)
pdf.output(f'{prot_pdb}.pdf')

# Eliminar archivos temporales de gráficos
os.remove(glycine_plot_name)
os.remove(proline_plot_name)
os.remove(plot_name)
os.remove(hist_plot_name)
os.remove(classification_name)
