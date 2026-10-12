#Proyecto: Filogenia de genotipos del Virus de Hepatitis B (NCBI real)

from Bio import Entrez, SeqIO

#Correo global para la API de NCBI
Entrez.email = "geovajfb@gmail.com"

#Datos reales de NCBI (números de acceso de genotipos de referencia):
accesiones = {
    "Genotipo_A": "X02763",
    "Genotipo_B": "D00329",
    "Genotipo_C": "AB033554",
    "Genotipo_D": "X65257",
    "Genotipo_E": "X75657"
}

#Desarrollo del proyecto:
#Descargar las secuencias desde NCBI

# 1. Usa Bio.Entrez para conectarte a NCBI y descargar cada secuencia por su código de acceso
def descargar_secuencia(accession):
    handle = Entrez.efetch(
        db="nucleotide", 
        id=accession, 
        rettype="fasta", 
        retmode="text"
    )

# 2. Leer y parsear el formato FASTA
    secuencia = SeqIO.read(handle, "fasta")
    handle.close()
    return secuencia

#Construir el diccionario de secuencias
secuencias_hbv ={genotipo: descargar_secuencia(acc) for genotipo, acc in accesiones.items()}

#Verificar la descarga
for genotipo, secuencia in secuencias_hbv.items():
    print(f"{genotipo} ({secuencia.id}): {len(secuencia.seq)} pb descargados")
    
#Crear una lista para guardar las secuencias con IDs limpios
secuencia_para_guardar =[]

# Asignar un ID limpio a cada registro (ej: "Genotipo_A")
for genotipo, secuencia in secuencias_hbv.items():
    secuencia.id = genotipo
    secuencia.description = f"HBV {genotipo}"
    secuencia_para_guardar.append(secuencia)
    
# Guardar en disco local como archivo multi-FASTA
SeqIO.write(secuencia_para_guardar, "secuencias_HBV.fasta", "fasta")
print("\n Archivo 'secuencias_HBV.fasta' guardado exitosamente.")

#ImportaR los módulos y define la función que tomará tu lista de secuencias:
from Bio.Align import MultipleSeqAlignment, PairwiseAligner


from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def alinear_secuencias(secuencias):
    #Configurar el alineador de pares
    alineador = PairwiseAligner()
    alineador.mode = "global"
    alineador.open_gap_score = -10
    alineador.extend_gap_score = -0.5
    
    # 2. La primera secuencia será la referencia base
    ref_secuencia = secuencias[0]
    secuencias_alineadas = [ref_secuencia] 
    
    # 3. Alinear las demás secuencias contra la referencia
    for secuencia in secuencias[1:]:
        alineamiento = alineador.align(ref_secuencia.seq, secuencia.seq)
        mejor_alineamiento = alineamiento[0]
        
        # Extraer la secuencia problema alineada (con sus gaps '-')
        cadena_alineada = mejor_alineamiento[1]
        
        secuencia_alineada = SeqRecord(
            Seq(cadena_alineada),
            id=secuencia.id,
            description=secuencia.description
        )
        secuencias_alineadas.append(secuencia_alineada)
        
        # 4. Ajustar longitudes para asegurar matriz rectangular perfecta
    max_len = max(len(s.seq) for s in secuencias_alineadas)
    for s in secuencias_alineadas:
        if len(s.seq) < max_len:
            s.seq = Seq(str(s.seq) + "-" * (max_len - len(s.seq)))
                
    return MultipleSeqAlignment(secuencias_alineadas)


# Llámar a la función pasando la lista de secuencias descargadas
alineamiento_hbv = alinear_secuencias(secuencia_para_guardar)

# Imprimir información del alineamiento
print("\n Alineamiento Multiple Completado!")
print(f"Numero de secuencias en el alimiento: {len(alineamiento_hbv)}")
print(f"Longitud alineada final (con gaps): {alineamiento_hbv.get_alignment_length()} pb")

#CONSTRUIR EL ARBOL FILOGENETICO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

# Crear el calculador de distancias
def construir_arbol_filogenetico(alineamiento):
    calculator = DistanceCalculator("identity")
    matriz_distancias = calculator.get_distance(alineamiento)

    #Crear el constructor de arboles filogeneticos
    constructor = DistanceTreeConstructor()

    #Contruir el arbol con Neighbor-Joining (NJ) y mostrarlo
    arbol = constructor.nj(matriz_distancias)
    
    # Organizar el árbol para una mejor visualización
    arbol.ladderize()
    return arbol  

## --- EJECUCION ---
arbol_filogenetico_hbv = construir_arbol_filogenetico(alineamiento_hbv)

#Mostrar el arbol filogenetico en consola
import matplotlib.pyplot as plt
from Bio import Phylo

print("\n--- Arbol Filogenetico de Genotipo de HBV (Neighbor-Joining) ---")
Phylo.draw_ascii(arbol_filogenetico_hbv)

# Guardar el arbol filogenetico en un archivo de imagen (PNG)
fig = plt.figure(figsize=(10, 6))
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(arbol_filogenetico_hbv, axes=axes, do_show=False)
plt.savefig("Arbol_Filogenetico_HBV.png", dpi=300)
print("\n Arbol filogenetico guardado como: 'Arbol_Filogenetico_HBV.png' exitosamente.")

#OPCIONAL: Mostrar el arbol filogenetico en una ventana gráfica
#plt.show()

#CONSTRUIR ARBOL FILOGENETICO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

def construir_arbol_filogenetico(alineamiento):
    calculator = DistanceCalculator("identity")
    matriz_distancias = calculator.get_distance(alineamiento)
    constructor = DistanceTreeConstructor()
    arbol = constructor.nj(matriz_distancias)
    arbol.ladderize()
    
    return arbol, matriz_distancias

#CREAR FUNCION PARA EL REPORTE FINAL DEL PROYECTO

import pandas as pd


def generar_reporte_final(accesiones):
    print("=" * 55)
    print("ANALISIS FILOGENÉTICO - GENOTIPOS DEL VIRUS DE HEPATITIS B (HBV)")
    print("=" * 55)
    
    #DESCARGAR SECUENCIAS DE NCBI
    print("\nDescargando sencuencias de NCBI...")
    secuencias_hbv = {genotipo: descargar_secuencia(acc) for genotipo, acc in accesiones.items()}
    
    for genotipo, secuencia in secuencias_hbv.items():
        print(f"{genotipo} ({secuencia.id}): {len(secuencia.seq)} pb descargados")
        
    #LIMPIAR IDs Y GUARDAR EN ARCHIVO FASTA
    secuencia_para_guardar = []
    for genotipo, secuencia in secuencias_hbv.items():
        secuencia.id = genotipo
        secuencia.description = f"HBV {genotipo}"
        secuencia_para_guardar.append(secuencia)
        
    SeqIO.write(secuencia_para_guardar, "secuencias_HBV.fasta", "fasta")
    print("\n Archivo 'secuencias_HBV.fasta' guardado exitosamente.")
    
    # PASO 3: ALINEAR SECUENCIAS
    print("\nAlienando secuencias...")
    alineamiento_hbv = alinear_secuencias(secuencia_para_guardar)
    print(f"Secuencias alineadas: {len(alineamiento_hbv)}")
    print(f"Longitud alineada (con gaps): {alineamiento_hbv.get_alignment_length()} pb")
    
    #PASO 4: CONSTRUIR ARBOL FILOGENETICO
    print("\n Construyendo arbol filogenetico (Neighbor - Joining)...")
    arbol_filogenetico_hbv, matriz_distancias = construir_arbol_filogenetico(alineamiento_hbv)
    
    print("\n Arbol Filogenetico (vista en consola) ---")
    Phylo.draw_ascii(arbol_filogenetico_hbv)
    
    #PASO 5: EXPORTAR MATRICES DE DISTANCIAS A CSV
    nombres = matriz_distancias.names
    df_distancias = pd.DataFrame(
        [[matriz_distancias[i, j] for j in nombres] for i in nombres],
        index=nombres,
        columns=nombres
    )
    df_distancias.to_csv("Matriz_Distancias_HBV.csv")
    print("\n Matriz de distancias guardada como 'Matriz_Distancias_HBV.csv' exitosamente.")
    
    
    # PASO 6: GUARDAR GRAFICO DEL ARBOL FILOGENETICO
    fig = plt.figure(figsize=(10,6))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(arbol_filogenetico_hbv, axes=axes, do_show=False)
    plt.savefig("Arbol_Filogenetico_HBV.png", dpi=300)
    print("\n Arbol filogenetico guardado como: ' Arbol_Filogenetico_HBV.png' exitosamente.")
    
    print("\n" + "=" * 55)
    print( "ANALISIS COMPLETADO - GENOTIPOS DEL VIRUS DE HEPATITIS B (HBV)")
    print("=" * 55)
    
    return{
        "secuencias": secuencias_hbv,
        "alineamiento": alineamiento_hbv,
        "arbol": arbol_filogenetico_hbv,
        "matriz_distancias": df_distancias
    }
    
#EJECUCION DE LA FUNCION PARA GENERAR EL REPORTE FINAL

reporte_final = generar_reporte_final(accesiones)
plt.show()
