import xmltodict as xtd
import folium
import numpy as np
import time
from IPython.display import display
from IPython.display import clear_output 
with open('mapHSR.osm', "rb") as osm_fn:
  map_osm = xtd.parse(osm_fn)['osm']
     
  map_osm.keys()


ymax = map_osm['bounds']['@maxlat']
ymin = map_osm['bounds']['@minlat']
xmax = map_osm['bounds']['@maxlon']
xmin = map_osm['bounds']['@minlon']

parsed_bounds = [xmin, xmax, ymin, ymax]



Node=map_osm['node']

Nnodes=len(Node)

Nodeid = [0]*Nnodes
xy = []

for i in range(Nnodes):
  Nodeid[i]=float(Node[i]['@id'])
  x=float(Node[i]['@lat'])
  y=float(Node[i]['@lon'])
  xy.append([x,y])

parsed_node={'id':Nodeid, 'xy':xy}

     



Way=map_osm['way']

Nways=len(Way)

Wayid=[0]*Nways
nodes_in_way=[0]*Nways
tags=[0]*Nways
for i in range(Nways):
  tempWay = Way[i]
  
  Wayid[i] = float(tempWay['@id'])

  Nnd=len(tempWay['nd'])
  ndTemp=[0]*Nnd
  
  for j in range(Nnd):
    ndTemp[j]=float(tempWay['nd'][j]['@ref'])
  
  nodes_in_way[i] = ndTemp

  if 'tag' in tempWay.keys():
    if type(tempWay['tag']) is list:
      tags[i]=tempWay['tag']
    else:
      tags[i]=[tempWay['tag']]
  else:
    tags[i]=[]

parsed_way={'id':Wayid,'nodes':nodes_in_way, 'tags':tags}



Relation=map_osm['relation']

Nrelation=len(Relation)

Relationid=[0]*Nrelation

for i in range(Nrelation):
  currentRelation = Relation[i]
  currentId=currentRelation['@id']
  Relationid[i]=float(currentId)

parsed_relation={'id':Relationid}



parsed_osm={
    'bounds':parsed_bounds,
    'way':parsed_way,
    'node':parsed_node,
    'relation':parsed_relation,
    'attributes':map_osm.keys()
}
     


bounds=parsed_osm['bounds']
way=parsed_osm['way']
node=parsed_osm['node']
relation=parsed_osm['relation']
     

ways_num = len(way['id'])
ways_node_set=way['nodes']
node_ids = dict()
n = len(node['id'])
for i in range(n):
  node_ids[node['id'][i]] = i

road_vals = ['highway', 'motorway', 'motorway_link', 'trunk', 'trunk_link',
             'primary', 'primary_link', 'secondary', 'secondary_link',
             'tertiary', 'road', 'residential', 'living_street',
             'service', 'services', 'motorway_junction']
     

def create_connectivity():
  connectivity_matrix = np.full((Nnodes,Nnodes), float('inf'))
  np.fill_diagonal(connectivity_matrix, 0)

  for currentWay in range(ways_num):
    skip = True
    for i in way['tags'][currentWay]:
      if i['@k'] in road_vals:
        skip = False
        break
    
    if skip:
      continue

    nodeset=ways_node_set[currentWay]
    nodes_num=len(nodeset)

    currentWayID = way['id'][currentWay]

    for firstnode_local_index in range(nodes_num):
      firstnode_id = nodeset[firstnode_local_index]
      firstnode_index = node_ids.get(firstnode_id, -1)
      if firstnode_index==-1: continue 

      for othernode_local_index in range(firstnode_local_index+1, nodes_num):
        othernode_id=nodeset[othernode_local_index]
        othernode_index = node_ids.get(othernode_id, -1)
        if othernode_index==-1: continue 

        if(firstnode_id != othernode_id and connectivity_matrix[firstnode_index,othernode_index]==float('inf')):
          connectivity_matrix[firstnode_index, othernode_index] = 1
          connectivity_matrix[othernode_index, firstnode_index] = 1

  return connectivity_matrix
     
def dijkstra(source, connectivity_matrix, p):
  s = dict()
  s[source] = True
  p[source] = source

  v = len(connectivity_matrix)
  u = source
  d_u = float('inf')
  for i in range(v):
    if i != source and connectivity_matrix[source][i] < d_u:
      u = i
      d_u = connectivity_matrix[source][i]
  s[u] = True
  p[u] = source

  i = v-2
  while i > 0:
    u_x = source
    d_u = float('inf')

    for j in range(v):
      if s.get(j, False) == False and connectivity_matrix[source][u] != float('inf') and connectivity_matrix[u][j] != float('inf'):
        k = connectivity_matrix[source][u] + connectivity_matrix[u][j]
        connectivity_matrix[source][j] = min(connectivity_matrix[source][j], k)
        connectivity_matrix[j][source] = connectivity_matrix[source][j]

        if connectivity_matrix[source][j] == k:
          p[j] = u
        elif connectivity_matrix[source][j] == 1:
          p[j] = source

        if connectivity_matrix[source][j] < d_u:
          u_x = j
          d_u = connectivity_matrix[source][j]

    if u_x == source:
      break

    s[u_x] = True
    u = u_x
    i -= 1
     


def plot_routes(s, connectivity_matrix):
  p = dict()
  dijkstra(s, connectivity_matrix, p)

  nodes_routes_values=[]
  for i in p.keys():
    adder=[i,0]
    while p[i] != i:
      adder[1]+=1
      i = p[i]
    nodes_routes_values.append(adder)

  return nodes_routes_values,p
     


def build_map(i,p):
  node_cds = [(node['xy'][i][0], node['xy'][i][1])]
  while p[i] != i:
    node_cds.append((node['xy'][p[i]][0], node['xy'][p[i]][1]))
    i = p[i]

  map = folium.Map(location = node_cds[-1], zoom_start = 15)

  folium.CircleMarker(node_cds[-1], radius=5, color="blue", fill=True, fill_color="orange").add_to(map)
  folium.Marker(node_cds[0], icon = folium.Icon(color="blue", icon="circle", prefix='fa')).add_to(map)
    
  folium.PolyLine(locations = node_cds, weight=5, color="blue", opacity="0.75", dash_array=10).add_to(map)
    
  return map
     


connectivity_matrix = create_connectivity()
nodes_routes_values,p = plot_routes(125, connectivity_matrix)
     

print(nodes_routes_values)
map = build_map(298,p)
display(map)

'''for i,j in nodes_routes_values:
  if j > 2:
    map = build_map(i,p)
    display(map)
    time.sleep(2)
    clear_output()'''