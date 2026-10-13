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

Cómo leer los resultados

1. B y C son los más parecidos entre sí (5.4 % de posiciones distintas) y en el árbol salen de una misma rama: forman un grupo.
2. D y E también forman un grupo (6.6 % de diferencias).
3. A no se agrupa con ninguno de los dos pares. Se diferencia de todos los demás entre 8.2 % y 9.7 %, casi por igual, así que el método no logra ubicarlo cerca de un grupo en particular.
4. El árbol es sin raíz: lo que importa es qué genotipos comparten una rama, no si quedan arriba o abajo en el dibujo.
5. Las diferencias van de ~5 % a ~11 %, el orden de magnitud esperado entre genotipos del mismo virus (la definición clásica habla de más de 8 % en el genoma completo). Mis valores para B–C y D–E quedan algo por debajo de ese umbral, probablemente porque uso una distancia simple (proporción de posiciones distintas, sin corrección por sustituciones múltiples) y un solo genoma de referencia por genotipo.

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

- `Proyecto_Virologia_HBV.py`: script completo
- `secuencias_HBV.fasta`: secuencias descargadas
- `Matriz_Distancias_HBV.csv`: matriz de distancias
- `Arbol_filogenetico_HBV.png`: árbol filogenético

## Qué aprendí

- Descargar y manejar datos públicos con Biopython (`Entrez`, `SeqIO`).
- Detectar un problema real de datos: al principio las distancias salían de ~75 %, lo cual es imposible entre genotipos del mismo virus. La causa era que el genoma A empezaba en otro punto del genoma circular. Validar que los resultados tengan sentido biológico fue tan importante como que el código corriera.
- Construir y visualizar árboles con `Bio.Phylo`.

## Limitaciones

- Con solo cinco secuencias y un único método (Neighbor-Joining sin bootstrap), el árbol no tiene soporte estadístico en sus ramas.
- Las distancias salen de alineamientos por pares, no de un alineamiento múltiple (MAFFT, Clustal). Para un análisis riguroso conviene un MSA y modelos de sustitución más completos (máxima verosimilitud o bayesiano).
- Solo se usa una secuencia de referencia por genotipo.
