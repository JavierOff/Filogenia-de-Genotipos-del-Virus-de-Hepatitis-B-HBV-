# Filogenia de Genotipos del Virus de Hepatitis B (HBV)

Pipeline en Python que descarga genomas reales de HBV desde NCBI, calcula las distancias genéticas entre cinco genotipos (A–E) y construye un árbol filogenético con el método Neighbor-Joining.

![Árbol filogenético de HBV](Arbol_Filogenetico_HBV.png)

## Sobre el proyecto

Como biólogo aprendiendo bioinformática, quise replicar un flujo de trabajo real de virología molecular: obtener secuencias públicas, compararlas y visualizar sus relaciones evolutivas. Elegí HBV porque su genoma es pequeño (~3.2 kb), y porque sus genotipos son un caso clásico de estudio en epidemiología molecular.

## Resultados

| | A | B | C | D | E |
|---|---|---|---|---|---|
| **A** | 0.000 | 0.082 | 0.091 | 0.094 | 0.097 |
| **B** | 0.082 | 0.000 | 0.054 | 0.097 | 0.101 |
| **C** | 0.091 | 0.054 | 0.000 | 0.098 | 0.107 |
| **D** | 0.094 | 0.097 | 0.098 | 0.000 | 0.066 |
| **E** | 0.097 | 0.101 | 0.107 | 0.066 | 0.000 |

*Proporción de posiciones diferentes entre cada par de genomas.*

- Los genotipos **B y C** forman un grupo (distancia 5.4 %), y **D y E** forman otro (6.6 %).
- El genotipo **A** queda separado de ambos grupos.
- Las distancias entre genotipos están entre ~5 % y ~11 %, un rango coherente con la divergencia descrita para los genotipos de HBV (>8 % entre genotipos distintos).

## Cómo funciona

1. **Descarga desde NCBI** con `Bio.Entrez` (accesiones reales de genotipos de referencia) y exportación a `secuencias_HBV.fasta`.
2. **Origen común**: el genoma de HBV es circular y los registros de NCBI no siempre empiezan en el mismo punto (el genotipo A estaba desplazado ~1400 pb). El script detecta el desplazamiento por k-meros compartidos y rota cada secuencia al mismo origen.
3. **Distancias**: alineamiento global por pares (`Bio.Align.PairwiseAligner`) y distancia-p entre cada par de secuencias.
4. **Árbol**: Neighbor-Joining con `Bio.Phylo.TreeConstruction`.
5. **Salidas**: matriz en CSV, árbol en consola (ASCII) y en PNG.

## Datos

| Genotipo | Accesión NCBI |
|----------|---------------|
| A | X02763 |
| B | D00329 |
| C | AB033554 |
| D | X65257 |
| E | X75657 |

## Uso

```bash
pip install biopython pandas matplotlib
python Proyecto_Virologia.py
```

Antes de ejecutar, cambia `Entrez.email` por tu propio correo (NCBI lo solicita para identificar las consultas).

## Archivos

- `Proyecto_Virologia.py`: script completo
- `secuencias_HBV.fasta`: secuencias descargadas
- `Matriz_Distancias_HBV.csv`: matriz de distancias
- `Arbol_Filogenetico_HBV.png`: árbol filogenético

## Qué aprendí

- Descargar y manejar datos públicos con Biopython (`Entrez`, `SeqIO`).
- Detectar un problema real de datos: al principio las distancias salían de ~75 %, lo cual es imposible entre genotipos del mismo virus. La causa era que el genoma A empezaba en otro punto del genoma circular. Validar que los resultados tengan sentido biológico fue tan importante como que el código corriera.
- Construir y visualizar árboles con `Bio.Phylo`.

## Limitaciones

- Con solo cinco secuencias y un único método (Neighbor-Joining sin bootstrap), el árbol no tiene soporte estadístico en sus ramas.
- Las distancias salen de alineamientos por pares, no de un alineamiento múltiple (MAFFT, Clustal). Para un análisis riguroso conviene un MSA y modelos de sustitución más completos (máxima verosimilitud o bayesiano).
- Solo se usa una secuencia de referencia por genotipo.
