# ToBaCCo

Aqui agregare alguna notas sobre el uso de tobacco para generar MOF.

## Instalacion

requerimientos:
```
numpy
networkx
scipy
```
luego se desgarga el repositorio
```
git clone https://github.com/tobacco-mofs/tobacco_3.0.git
```
La ultima modificacion fue hecha hace 3 anos.

## Usarlo como modulo

Se ha creado una rama de git llamada testing para modificar y probar parte del codigo, tambien se ha agregado un archivo `__init__.py` en el repertorio base, tabmien se agrego al path de sitios de python:
```
echo "/home/ovillegas/gitproyects/tobacco_3.0" > /home/ovillegas/.local/lib/python3.12/site-packages/tobacco.pth
```
, esto hace que python lo reconozca como un modulo.

## Implementacion

El script principal es `tobacco.py`.


1) Este comienza con la inicializacion de diferentes configuraciones.`import configuration`, aqui abajo podemos ver las configuraciones que se crean por defecto:

```
IGNORE_ALL_ERRORS = False
PRINT = False
CONNECTION_SITE_BOND_LENGTH = 1.54
WRITE_CHECK_FILES = False
WRITE_CIF = True
ALL_NODE_COMBINATIONS = False
USER_SPECIFIED_NODE_ASSIGNMENT = False
COMBINATORIAL_EDGE_ASSIGNMENT = True
CHARGES = True
SYMMETRY_TOL = {2:0.10, 3:0.12, 4:0.35, 5:0.25, 6:0.45, 7:0.35, 8:0.40, 9:0.60, 10:0.60, 12:0.60}
BOND_TOL = 5.0
ORIENTATION_DEPENDENT_NODES = False
PLACE_EDGES_BETWEEN_CONNECTION_POINTS = True
RECORD_CALLBACK = False
OUTPUT_SCALING_DATA = True
FIX_UC = (0,0,0,0,0,0)
MIN_CELL_LENGTH = 5.0
OPT_METHOD = 'L-BFGS-B'
PRE_SCALE = 1.00
SCALING_ITERATIONS = 1
SINGLE_METAL_MOFS_ONLY = True
MOFS_ONLY = True
MERGE_CATENATED_NETS = True
RUN_PARALLEL = False
REMOVE_DUMMY_ATOMS = True
```

Podemos apreciar en esta opcion que es posible correr el programa en paralelo.

2) Al correr el codigo como script principal este va recorrer tres carpetas principales: `templates`, `nodes`, `edges`. 

- `templates` contiene la informacion de la topologia de un cristal en formato cif.

- `nodes` contiene la definicion de los nodos como un archivo cif.

- Igualmente `edges` contiene efinicion de edges como archivos cif.

3) Antes de comenzar `ToBaCCo` leera y ordenara por orden alfabetico todas las topologias presentes en el repertorio `templates`:

```
templates = sorted(os.listdir('templates'))
```

4) Luego dependiendo de si usara la opcion en paralelo o no ejecutara la funcion central:

```
if RUN_PARALLEL:
    run_tobacco_parallel(templates, CHARGES)
else:
    run_tobacco_serial(templates, CHARGES)
```

### La funcion `run_tobacco_serial`

Al correr esta funcion esta recibe como parametro la lista de topologias y un boleano indicando si se agregaran las cargas o no, por defecto es verdadero.

```
run_tobacco_serial(templates, CHARGES)
```

La funcion esta definida en las siguientes lineas:

```
def run_tobacco_serial(templates, CHARGES):

    if IGNORE_ALL_ERRORS:
        for template in templates:
            try:
                run_template(template)
            except Exception as e:
                print()
                print('*****************************************************************')
                print('ERROR for template :',template)      
                print('error message:',e)
                print('continuing to next template...')                          
                print('*****************************************************************')
                print()
    else:
        for template in templates:
            run_template(template)
```

Aqui vemos que para una topologia particular este ejecutara la funcion `run_template(template)`. Lo primero que hace esta funcion es tranformar la informacion de la topologia como un grapho y va iterar sobre los nodos que lo componen `for net in ct2g(template):`

### Correr ToBaCCo en paraelo

La funcion que gobierna este proceso es:

```
def run_tobacco_parallel(templates, CHARGES):
    
    print('running parallel on', multiprocessing.cpu_count(), 'processors...')

    args = [template  for template in templates]
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map_async(run_template, args) 
    pool.close()
    pool.join()
```

**NOTA: He observado que la topologia cds definida en Tobacco es diferente a la que define en RCSR**

Se correra un test usando todas las 2600 topologia del repertorio `template_database` y usando la funcion en paralelo.

### La funcion principal de Tobacco

La funcion principal de Tobacco es `run_template`, esta sera la generadora de MOF directamente a partir de una topologia.

La primera topologia tomada por tobacco es `aab`, donde se observan 19 nodos y tiene el grupo Pmmm (P=primitive, there is only 1 reticular point inside the cell, esto quiere decir que cada punto de la red es unico y no hay puntos en la cercania que se repitan, mmm=indica la simetria de planos de reflexion perpendiculares entre si)

Las topologia son registradas en formato `nx.MultiGraph()`, la funcion `ct2g` identificara diferentes lineas que comienzan por `_`

Se agregan los nodos y los edges, los nodos tendran la siguientes propiedades:
```
G.add_node(s[0], type=ty, index=nc, ccoords=c_nvec, fcoords=f_nvec, cn=[], cifname=[])
```
los edges son agregados como edges de grafos: 
```
G.add_edge(s[0],s[1], key=(ne,lbl[0],lbl[1],lbl[2]), label=lbl , length=le, fcoords=ef_coords, ccoords=ec_coords, index=ne, pd=(s[0],s[1]))
```

Esta funcion `ct2g` regresa un generador python, no regresa la estructura de grafo.
```
yield (SG, start, unit_cell, set(cns), set(e_types), cifname, aL, bL, cL, alpha, beta, gamma, max_le, catenation)
```

Regresando a la funcion `run_template`, se iterara sobre los nets presentes en una topologia dada:
```
for net in ct2g(template):
    TG, start, unit_cell, TVT, TET, TNAME, a, b, c, ang_alpha, ang_beta, ang_gamma, max_le, catenation = net
```


## Crear link simbolico
