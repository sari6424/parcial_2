# Análisis del Gráfico de Ramachandran de Proteínas

Un gráfico de Ramachandran es una representación gráfica de los ángulos de torsión phi (φ) y psi (ψ) de los residuos en una cadena polipeptídica. Estos ángulos permiten entender la estrucutra tridimensional de la  proteína, además describen la conformación de la cadena principal de la proteína. 

Este proyecto analiza los ángulos de torsión (phi y psi) de los residuos de proteínas a partir de un archivo PDB, genera el gráfico de Ramachandran general, uno de glicinas y otro de prolinas. También clasifica los residuos según su R y contabiliza los aminoácidos mostrados en el Ramachandran mediante un histograma.

## Características

- **Calcular Ángulos de Torsión**: Calcula los ángulos phi y psi para cada residuo en la proteína.
- **Generar Gráficos de Ramachandran**: Crea gráficos de Ramachandran para todos los residuos (a excepción del primero y el ultimo), residuos de prolina y residuos de glicina.
- **Clasificar Residuos**: Clasifica los residuos según su R: polar, polar positivo, polar negativo, no polar alifático y no polar aromático.
- **Generar Histogramas y Gráficos de Clasificación**: Visualiza la frecuencia de los residuos y sus clasificaciones.
- **Generar Informe en PDF**: Compila todos los gráficos generados en un informe en PDF.

## Instalación

Para ejecutar este proyecto se necesitan las siguientes bibliotecas:

- `numpy`
- `pandas`
- `seaborn`
- `matplotlib`
- `fpdf`
