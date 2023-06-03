###This is a source code with someome can create his own model with the 20_8_BrickUP element starting from GMSH modeler#####
import os
import openseespy.opensees as ops
import gmsh2opensees as g2o
import numpy as np
ops.wipe()
import gmsh
import math as mm
import time as tt

gmsh.initialize() #inizializzazione di gmesh
gmsh.open("Cubo_new.msh")

########################################################################################################################
####################### PVD ############################################################################################
########################################################################################################################


filename = 'cubottoplotto'
if not os.path.exists(filename):
    os.makedirs(filename)
ops.recorder('PVD', filename, 'disp', 'gausspoint', 'stress')
#ops.recorder('MPC', filename,'stress')

def floatingNodes():
    connectedNodes = []
    for ele in ops.getEleTags():
        for nd in ops.eleNodes(ele):
            connectedNodes.append(nd)

    definedNodes = ops.getNodeTags()

    # Use XOR operator, ^
    return list(set(connectedNodes) ^ set(definedNodes))

PhysGr=g2o.get_physical_groups_map(gmsh.model)

ops.model("basicBuilder","-ndm",3,"-ndf",4)  

matTag =  1
E = 210e9 #Pa
nd = 3 #dimension
nu = 0.3
rho = 7300. # kg/m3
friction = 31.0         #friction angle
phaseTransform = 26.0   #phase transformation angle
G1 = 9.e4
B1 = 22.e4
gamma =    0.600      # Newmark integration parameter

dT =   0.01           # time step for analysis, does not have to be the same as accDt.
numSteps= 2500       # number of time steps
rhoS  =1.80          # saturated mass density
rhoF  =1.00          # fluid mass density

Bfluid =2.2e6        # fluid shear modulus
fluid1 =1            # fluid material tag
solid1 =10           # solid material tag
perm   =1.e-5    #permeability (m/s)
accGravity =9.81  #acceleration of gravity
perm1   =[perm/accGravity/rhoF]    # actual value used in computation
accMul = 2                    # acceleration multiplier
pi = 3.1415926535
inclination = 0
massProportionalDamping   =0.0
InitStiffnessProportionalDamping =0.003
gravityX =[accGravity*np.sin(inclination/180.0*pi)] # gravity acceleration in X direction
gravityY =0.0                                        # gravity acceleration in Y direction
gravityZ =[accGravity*np.cos(inclination/180.0*pi)]  # gravity acceleration in Z direction

ops.nDMaterial('PressureDependMultiYield02', 1, 3, rhoS, G1, B1, friction, 2.5, 80, 0.5, phaseTransform, 0.067, 0.23, 0.06, 0.27)

g2o.get_elements_and_nodes_in_physical_group("Solid", gmsh.model)
elementTags, nodeTags, elementName, elementNnodes = g2o.get_elements_and_nodes_in_physical_group("Solid", gmsh.model)

LNT=len(nodeTags)
nodeNum = 0
nodeTot = []
nodeVec=[] #vettore corrispondente agli elementi
mem=[]
nodes_dic={}
x=0
y=0
z=0
for i in range(0,LNT):
    for j in range(0,20):
        nodeNum = nodeTags[i][j]
        nodeTot.append(nodeTags[i][j])
        x = gmsh.model.mesh.getNode(nodeNum)[0][0]
        y = gmsh.model.mesh.getNode(nodeNum)[0][1]
        z = gmsh.model.mesh.getNode(nodeNum)[0][2]
        nodes_dic[nodeNum] = (nodeNum, x, y, z)

nodeNum = 0

nodeEl = []
nodeVec = []  
mem = []

x = 0
y = 0
z = 0
for i in range(0,LNT):
    for j in range(0,8):
        nodeNum = nodeTags[i][j]
        nodeEl.append(nodeNum)
        mem.append(nodeNum)
        if j == 7:
            nodeVec.append(mem)
            mem = []
            j = 0

nodeVecU = np.unique(np.array(nodeVec, dtype=int).reshape(-1))
np.savetxt('Noderecord3D.txt',nodeVecU)
nodeNum = np.array(nodeEl)
nodeNumU= np.unique(np.array(nodeEl, dtype=int).reshape(-1))
np.savetxt('Noderecord4D.txt',nodeNumU)
xCoord = 0
yCoord = 0
zCoord = 0
Coord=[]
k=0
nodeNum=0

for nodeNum in nodeNumU:
    xCoord =  gmsh.model.mesh.getNode(nodeNum)[0][0]
    yCoord =  gmsh.model.mesh.getNode(nodeNum)[0][1]
    zCoord =  gmsh.model.mesh.getNode(nodeNum)[0][2]
    Coord_l=[round(xCoord,2),round(yCoord,2),round(zCoord,2)]
    ops.node(int(nodeNum), *Coord_l)                         
    

nodeTagsA = np.array(nodeTot)
nodeTagsU=np.unique(nodeTagsA)

np.savetxt('nodesInfo.txt',nodeTagsU)
xCoord = 0
yCoord = 0
zCoord = 0
Coord=[]
k=0
n = open("nodesInfo.dat","w")
for i in nodeTagsU:
    xCoord =  gmsh.model.mesh.getNode(i)[0][0]
    yCoord =  gmsh.model.mesh.getNode(i)[0][1]
    zCoord =  gmsh.model.mesh.getNode(i)[0][2]
    n.write(f"{i} {xCoord}    {yCoord}    {zCoord}\n")
n.close()

g2o.get_physical_groups_map(gmsh.model) 
TagPh = PhysGr['Top'][1]                                            
DimPh = PhysGr['Top'][0]
entities_top=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh)  
surfTags = gmsh.model.mesh.getElements(DimPh, entities_top[0])  


k = 0
surfEl = []
for k in surfTags[1][0]:
    surfEl.append(k)

surfNodeTags=[]
surfNodeCorn = []
surfNodeEdge = []
NodeListTop=[]
mem=[]
n = 0
i = 0
j = 0

len1 = len(surfTags[2][0])
for i in range(0, len1):
        n = surfTags[2][0][i]
        surfNodeTags.append(n)
        NodeListTop.append(n)
        mem.append(n)
        j=j+1
        if j ==4:
            surfNodeCorn.append(mem)
            mem=[]
        elif j == 8:
            surfNodeEdge.append(mem)
            j = 0
            mem=[]


xminSurf,yminSurf,zminSurf,xmaxSurf,ymaxSurf,zmaxSurf=gmsh.model.get_bounding_box(2,entities_top[0]) 

surfNodeCornTags=[]
LLA=len(surfNodeCorn)
LLB=len(surfNodeCorn[0])
for i in range(0,LLA):
    for j in range(0,LLB):
        surfNodeCornTags.append(surfNodeCorn[i][j])


surfNodeCornA = np.array(surfNodeCornTags,dtype=int).reshape(-1)
surfNodeCornU = np.unique(surfNodeCornA)
xCoordCorn = 0
yCoordCorn = 0
zCoordCorn = 0
Coord_Corn_l = []
CoordCorn={}
len2=0
nodi_bordo1=[]
nodi_bordo2=[]
nodi_bordo3=[]
nodi_bordo4=[]
i = 0
j = 0

for j in surfNodeCornU:
        xCoordCorn = gmsh.model.mesh.getNode(j)[0][0]
        yCoordCorn = gmsh.model.mesh.getNode(j)[0][1]
        zCoordCorn = gmsh.model.mesh.getNode(j)[0][2]
        Coord_Corn_l = [xCoordCorn, yCoordCorn, zCoordCorn]
        CoordCorn[j] = Coord_Corn_l
        if (xCoordCorn == xmaxSurf) or (xCoordCorn == xmaxSurf and (yCoordCorn == ymaxSurf or yCoordCorn == yminSurf )):
           ops.fix(int(j),1,0,0,1)
           nodi_bordo1.append(j)
           
        elif (xCoordCorn == xminSurf) or (xCoordCorn == xminSurf and (yCoordCorn == ymaxSurf or yCoordCorn == yminSurf )):
           ops.fix(int(j),1,0,0,1)
           nodi_bordo2.append(j)
       
        elif (yCoordCorn == ymaxSurf) or (yCoordCorn == ymaxSurf and (xCoordCorn == xmaxSurf or xCoordCorn == xminSurf )):
           ops.fix(int(j),0,1,0,1)
           nodi_bordo3.append(j)
           
        elif (yCoordCorn == yminSurf) or (yCoordCorn == yminSurf and (xCoordCorn == xmaxSurf or xCoordCorn == xminSurf )):
           ops.fix(int(j),0,1,0,1)
           nodi_bordo4.append(j)
           
        else:
           ops.fix(int(j),0,0,0,1)

g2o.get_physical_groups_map(gmsh.model)  
TagPh = PhysGr['Down'][1]
DimPh = PhysGr['Down'][0]
entities_down=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh)
surfTags = gmsh.model.mesh.getElements(DimPh, entities_down[0])  
entTagDown,elTagsDown,nodeTagsDown=gmsh.model.mesh.getElements(2, entities_down[0])

surfNodeTags=[]
surfNodeCorn = []
surfNodeEdge = []
NodeListDown = []
mem=[]
n = 0
i = 0
j = 0

len1 = len(surfTags[2][0])
for i in range(0, len1):
        n = surfTags[2][0][i]
        surfNodeTags.append(n)
        NodeListDown.append(n)
        mem.append(n)
        j=j+1
        if j ==4:
            surfNodeCorn.append(mem)
            mem=[]
        elif j == 8:
            surfNodeEdge.append(mem)
            j = 0
            mem=[]

xminSurf,yminSurf,zminSurf,xmaxSurf,ymaxSurf,zmaxSurf=gmsh.model.get_bounding_box(2,entities_down[0])

surfNodeCornTags=[]
LLA=len(surfNodeCorn)
LLB=len(surfNodeCorn[0])
for i in range(0,LLA):
    for j in range(0,LLB):
        surfNodeCornTags.append(surfNodeCorn[i][j])

surfNodeCornA = np.array(surfNodeCornTags,dtype=int).reshape(-1)
surfNodeCornU = np.unique(surfNodeCornA)
xCoordCorn = 0
yCoordCorn = 0
zCoordCorn = 0
Coord_Corn_l = []
CoordCorn={}
len2=0
nodi_bordo1=[]
nodi_bordo2=[]
nodi_bordo3=[]
nodi_bordo4=[]
i = 0
j = 0

for j in surfNodeCornU:
        tcl = open('cubotto.tcl', 'a')
        xCoordCorn = gmsh.model.mesh.getNode(j)[0][0]
        yCoordCorn = gmsh.model.mesh.getNode(j)[0][1]
        zCoordCorn = gmsh.model.mesh.getNode(j)[0][2]
        Coord_Corn_l = [xCoordCorn, yCoordCorn, zCoordCorn]
        CoordCorn[j] = Coord_Corn_l
        if (xCoordCorn == xmaxSurf) or (xCoordCorn == xmaxSurf and (yCoordCorn == ymaxSurf or yCoordCorn == yminSurf )):
           ops.fix(int(j),1,0,1,0)
           nodi_bordo1.append(j)

        elif (xCoordCorn == xminSurf) or (xCoordCorn == xminSurf and (yCoordCorn == ymaxSurf or yCoordCorn == yminSurf )):
           ops.fix(int(j),1,0,1,0)
           nodi_bordo2.append(j)

        elif (yCoordCorn == ymaxSurf) or (yCoordCorn == ymaxSurf and (xCoordCorn == xmaxSurf or xCoordCorn == xminSurf )):
           ops.fix(int(j),0,1,1,0)
           nodi_bordo3.append(j)

        elif (yCoordCorn == yminSurf) or (yCoordCorn == yminSurf and (xCoordCorn == xmaxSurf or xCoordCorn == xminSurf )):
           ops.fix(int(j),0,1,1,0)
           nodi_bordo4.append(j)

        else:
           ops.fix(int(j),0,0,1,0)

g2o.get_physical_groups_map(gmsh.model)  
TagPh = PhysGr['Solid'][1]
DimPh = PhysGr['Solid'][0]
entities_top=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh) 
entVol,VolTag,nodeTags = gmsh.model.mesh.getElements(DimPh, entities_top[0]) 

ops.model("basicBuilder","-ndm",3,"-ndf",3)                                 

j = 0
mem=[]
nodeTagsVec=[]
for i in nodeTags[0]:
    mem.append(i)
    j = j +1
    if j ==20:
        nodeTagsVec.append(mem)
        mem=[]
        j = 0

LNT=len(nodeTagsVec)
nodeNum2 = 0
nodeEl2=[]
nodeVec2=[]  
mem=[]
for i in range(0,LNT):
    for j in range(8,20):
        nodeNum2 = nodeTagsVec[i][j]
        nodeEl2.append(nodeNum2)
        mem.append(nodeNum2)
        if j == 19:
            nodeVec2.append(mem)
            mem = []
            j = 0

nodeNum2 = np.array(nodeEl2)

nodeNumU2= np.unique(np.array(nodeEl2, dtype=int).reshape(-1))
xCoord2 = 0
yCoord2 = 0
zCoord2 = 0
Coord2=[]
Coord_l2=[]
k=0

for nodeNum2 in nodeNumU2:
    xCoord2 =  gmsh.model.mesh.getNode(nodeNum2)[0][0]
    yCoord2 =  gmsh.model.mesh.getNode(nodeNum2)[0][1]
    zCoord2 =  gmsh.model.mesh.getNode(nodeNum2)[0][2]
    Coord_l2=[round(xCoord2,2),round(yCoord2,2),round(zCoord2,2)]
    ops.node(int(nodeNum2), *Coord_l2)                       

g2o.get_physical_groups_map(gmsh.model) 
TagPh = PhysGr['FixX'][1]
DimPh = PhysGr['FixX'][0]
entities_x=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh)
surfTags = gmsh.model.mesh.getElements(DimPh, entities_x[0])  
entTagDown,elTagsDown,nodeTagsDown=gmsh.model.mesh.getElements(2, entities_x[0])

xminSurf,yminSurf,zminSurf,xmaxSurf,ymaxSurf,zmaxSurf=gmsh.model.get_bounding_box(2,entities_x[0])

surfNodeTags=[]
surfNodeCorn = []
surfNodeEdge = []
NodeListX = []
mem=[]
n = 0
i = 0
j = 0

len1 = len(surfTags[2][0])
for i in range(0, len1):
        n = surfTags[2][0][i]
        surfNodeTags.append(n)
        NodeListX.append(n)
        mem.append(n)
        j=j+1
        if j == 4:
            surfNodeCorn.append(mem)
            mem=[]
        elif j == 8:
            surfNodeEdge.append(mem)
            j = 0
            mem=[]

i=0
LED=len(surfNodeEdge)
surfNodeCornA =[]
surfNodeEdgeA = np.array(surfNodeEdge)
surfNodeEdgeU = np.unique(surfNodeEdgeA)

LED=len(surfNodeEdge)
vecEdge=[]
for j in range(0,LED):
    for i in surfNodeEdge[j]:
        vecEdge.append(int(i))

vecEdgeA = np.array(vecEdge, dtype=int).reshape(-1)
vecEdgeU = np.unique(vecEdgeA)

for i in vecEdgeU:
    ops.fix(int(i),1,0,0)


g2o.get_physical_groups_map(gmsh.model)  
TagPh = PhysGr['FixY'][1]
DimPh = PhysGr['FixY'][0]
entities_y=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh)
surfTags = gmsh.model.mesh.getElements(DimPh, entities_y[0])  
entTagDown,elTagsDown,nodeTagsDown=gmsh.model.mesh.getElements(2, entities_y[0])

NodeListY=[]
surfNodeTags=[]
surfNodeCorn = []
surfNodeEdge = []
mem=[]
n = 0
i = 0
j = 0

len1 = len(surfTags[2][0])
for i in range(0, len1):
        n = surfTags[2][0][i]
        surfNodeTags.append(n)
        NodeListY.append(n)
        mem.append(n)
        j=j+1
        if j == 4:
            surfNodeCorn.append(mem)
            mem=[]
        elif j == 8:
            surfNodeEdge.append(mem)
            j = 0
            mem=[]


i=0
LED=len(surfNodeEdge)
surfNodeCornA =[]
surfNodeEdgeA = np.array(surfNodeEdge)
surfNodeEdgeU = np.unique(surfNodeEdgeA)

LED=len(surfNodeEdge)
vecEdge=[]
for j in range(0,LED):
    for i in surfNodeEdge[j]:
        vecEdge.append(int(i))

vecEdgeA = np.array(vecEdge, dtype=int).reshape(-1)
vecEdgeU = np.unique(vecEdgeA)

for i in vecEdgeU:
    ops.fix(int(i),0,1,0)
    
g2o.get_physical_groups_map(gmsh.model)  
TagPh = PhysGr['Down'][1]
DimPh = PhysGr['Down'][0]
entities_down=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh)
surfTags = gmsh.model.mesh.getElements(DimPh, entities_down[0]) 
entTagDown,elTagsDown,nodeTagsDown=gmsh.model.mesh.getElements(2, entities_down[0])

xminSurf,yminSurf,zminSurf,xmaxSurf,ymaxSurf,zmaxSurf=gmsh.model.get_bounding_box(2,entities_down[0])

surfNodeTags=[]
surfNodeCorn = []
surfNodeEdge = []
mem=[]
n = 0
i = 0
j = 0

len1 = len(surfTags[2][0])
for i in range(0, len1):
        n = surfTags[2][0][i]
        surfNodeTags.append(n)
        mem.append(n)
        j=j+1
        if j == 4:
            surfNodeCorn.append(mem)
            mem=[]
        elif j == 8:
            surfNodeEdge.append(mem)
            j = 0
            mem=[]


i=0
LED=len(surfNodeEdge)
surfNodeCornA =[]
surfNodeEdgeA = np.array(surfNodeEdge)
surfNodeEdgeU = np.unique(surfNodeEdgeA)

LED=len(surfNodeEdge)
vecEdge=[]
for j in range(0,LED):
    for i in surfNodeEdge[j]:
        vecEdge.append(int(i))

vecEdgeA = np.array(vecEdge, dtype=int).reshape(-1)
vecEdgeU = np.unique(vecEdgeA)

xCoord = 0
yCoord = 0
zCoord = 0
Coord_l=[]
Coord={}
for i in vecEdgeU:
    xCoord = gmsh.model.mesh.getNode(i)[0][0]
    yCoord = gmsh.model.mesh.getNode(i)[0][1]
    zCoord = gmsh.model.mesh.getNode(i)[0][2]
    Coord_l = [xCoord, yCoord, zCoord]
    Coord[i] = Coord_l
    if (xCoord == xmaxSurf) or (xCoord == xmaxSurf and (yCoord == ymaxSurf or yCoord == yminSurf)):
        nodi_bordo1.append(i)
        ops.fix(int(i),0,0,1)
                
    elif (xCoord == xminSurf) or (xCoord == xminSurf and (yCoord == ymaxSurf or yCoord == yminSurf)):
        nodi_bordo2.append(i)
        ops.fix(int(i), 0, 0, 1)
     

    elif (yCoord == ymaxSurf) or (yCoord == ymaxSurf and (xCoord == xmaxSurf or xCoord == xminSurf)):
        nodi_bordo3.append(i)
        ops.fix(int(i), 0, 0, 1)

       
    elif (yCoord == yminSurf) or (yCoord == yminSurf and (xCoord == xmaxSurf or xCoord == xminSurf)):
        nodi_bordo4.append(i)
        ops.fix(int(i), 0, 0, 1)

    else:
        ops.fix(int(i), 0, 0, 1)                                
       
ops.model("basicBuilder","-ndm",3,"-ndf",4)  

g2o.get_physical_groups_map(gmsh.model)  
TagPh = PhysGr['Solid'][1]
DimPh = PhysGr['Solid'][0]
entities_Vol=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh)
VolTags = gmsh.model.mesh.getElements(DimPh, entities_Vol[0])  
entTagVol,elTagsVol,nodeTagsVol=gmsh.model.mesh.getElements(DimPh, entities_Vol[0])

f=open("cubottoGID.msh","w")
el=open("elementInfo.dat","w")
f.write("MESH dimension 3 ElemType Hexahedra Nnode 20\n")  
f.write("Coordinates\n")
f.write("#node_number   coord_x   coord_y   coord_z\n")


for j in nodeTot:
    f.write(f"{j}   {gmsh.model.mesh.getNode(j)[0][0]}   {gmsh.model.mesh.getNode(j)[0][1]}   {gmsh.model.mesh.getNode(j)[0][2]}\n")
f.write("end coordinates\n")
f.write("Elements\n")

f.write("# element   nodo1  nodo2   nodo3   nodo4   nodo5   nodo6   nodo7   nodo8   nodo9   nodo10   nodo11   nodo12   nodo13   nodo14   nodo15   nodo16   nodo17   nodo18   nodo19   nodo20\n")
 

eleTag=[]
eleNodes = []
processed=[]
mem=[]
i = 0
j = 0
k = 0
LET = len(elementTags)
for i in range(0,LET):
    eleTag.append(elTagsVol[0][i]) 

for k in nodeTagsVol[0]:
    mem.append(k)
    j = j + 1
    if j == 20:
        eleNodes.append(mem)
        mem = []
        j = 0

tcl.close()
nodes=[]
connList={}
for i in range(0,LET):
    elem = int(eleTag[i])
    nodo1 = int(eleNodes[i][4])
    # print(nodes_dic[nodo1])
    nodo2 = int(eleNodes[i][5])
    # print(nodes_dic[nodo2])
    nodo3 = int(eleNodes[i][1])
    # print(nodes_dic[nodo3])
    nodo4 = int(eleNodes[i][0])
    # print(nodes_dic[nodo4])
    nodo5 = int(eleNodes[i][7])
    # print(nodes_dic[nodo5])
    nodo6 = int(eleNodes[i][6])
    # print(nodes_dic[nodo6])
    nodo7 = int(eleNodes[i][2])
    # print(nodes_dic[nodo7])
    nodo8 = int(eleNodes[i][3])
    # print(nodes_dic[nodo8])
    nodo9 = int(eleNodes[i][16])
    # print(nodes_dic[nodo9])
    nodo10 = int(eleNodes[i][12])
    # print(nodes_dic[nodo10])
    nodo11 = int(eleNodes[i][8])
    # print(nodes_dic[nodo11])
    nodo12 = int(eleNodes[i][10])
    # print(nodes_dic[nodo12])
    nodo13 = int(eleNodes[i][19])
    # print(nodes_dic[nodo13])
    nodo14 = int(eleNodes[i][14])
    # print(nodes_dic[nodo14])
    nodo15 = int(eleNodes[i][13])
    # print(nodes_dic[nodo15])
    nodo16 = int(eleNodes[i][15])
    # print(nodes_dic[nodo16])
    nodo17 = int(eleNodes[i][17])
    # print(nodes_dic[nodo17])
    nodo18 = int(eleNodes[i][18])
    # print(nodes_dic[nodo18])
    nodo19 = int(eleNodes[i][11])
    # print(nodes_dic[nodo19])
    nodo20 = int(eleNodes[i][9])
    # print(nodes_dic[nodo20])
    nodes_l = [nodo1,nodo2,nodo3,nodo4,nodo5,nodo6,nodo7,nodo8,nodo9, nodo10, nodo11,nodo12, nodo13, nodo14,nodo15,nodo16, nodo17,nodo18, nodo19, nodo20]
    connList[elem]=(nodes_l)
    f.write(f"{elem} {nodo1} {nodo2} {nodo3} {nodo4} {nodo5} {nodo6} {nodo7} {nodo8} {nodo9} {nodo10} {nodo11} {nodo12} {nodo17} {nodo18} {nodo19} {nodo20} {nodo13} {nodo14} {nodo15} {nodo16} \n")
    el.write(f"{elem} {nodo1} {nodo2} {nodo3} {nodo4} {nodo5} {nodo6} {nodo7} {nodo8} {nodo9} {nodo10} {nodo11} {nodo12} {nodo17} {nodo18} {nodo19} {nodo20} {nodo13} {nodo14} {nodo15} {nodo16}\n")
     
    ops.element('20_8_BrickUP',elem,*nodes_l, 1, 2.2e6, 1, 1.0, 1.0, 1.0, 0.0, 0.0,-20)

    
f.write("end elements")
f.close()
el.close()


ops.model("basicBuilder","-ndm",3,"-ndf",3)



dashX = 757

dashF = 758
ops.node(dashX,5.0,5.0,0.0)

ops.node(dashF,5.0,5.0,0.0)

ops.fix(dashF,1,1,1)

ops.fix(dashX,0,1,1)

NodeListDownAR=np.array(NodeListDown)
NodeListDownU = np.unique(NodeListDownAR)
SupDownDic={}
for i in NodeListDownU:
    xCoord = gmsh.model.mesh.getNode(i)[0][0]
    yCoord = gmsh.model.mesh.getNode(i)[0][1]
    zCoord = gmsh.model.mesh.getNode(i)[0][2]
    SupDownDic[i]=(xCoord,yCoord,zCoord)

ops.equalDOF(458, dashX,1)

colArea =   100 # l'elemento cubico ha lato 2
rockVS =       700.0
rockDen  =     2.5
dashpotCoeff = rockVS*rockDen*colArea
ops.uniaxialMaterial('Viscous',2,dashpotCoeff, 1)
Len=len(eleTag)
nElemT =eleTag[Len-1]

ops.element('zeroLength', int(nElemT+1), int(dashF), int(dashX), '-mat', 2, '-dir', 1)

g2o.get_physical_groups_map(gmsh.model) #dictionary - numero tag e dimensione ph ent.
TagPh = PhysGr['Solid'][1]
DimPh = PhysGr['Solid'][0]
entities_top=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh) # riportare dim e tag della superficie i cui i punti sono di interesse
entVol,VolTag,nodeTags = gmsh.model.mesh.getElements(DimPh, entities_top[0]) # punti della ph ent selezionata
nElemT = ops.getNumElements()

load_nodeList=np.loadtxt('nodesInfo.txt')
nodeList=[]

for i in range(len(load_nodeList)):
    nodeList.append(int(load_nodeList[i]))

load_nodeList4D=np.loadtxt('Noderecord4D.txt')
nodeList4D=[]

for i in range(len(load_nodeList4D)):
    nodeList4D.append(int(load_nodeList4D[i]))

ops.recorder('Node','-file','Gdisplacement.out','-time','-node',*nodeList,'-dof', 1, 2, 3, 'disp')
ops.recorder('Node','-file','Gacceleration.out','-time','-node',*nodeList,'-dof', 1, 2, 3, 'accel')
ops.recorder('Node','-file','GporePressure.out','-time','-node',*nodeList,'-dof', 4, 'vel')

ops.recorder('Element','-file','Gstress.out','-time','-eleRange', 211,335,'material','1','stress')
ops.recorder('Element','-file','Ggauss.out','-time','-eleRange', 211,335,'material','1','gausspoint')
ops.recorder('Element','-file','Gstrain.out','-time','-eleRange', 211,335,'material','1','strain')

motionDT = 0.005

motionSteps = 2500


#---RAYLEIGH DAMPING PARAMETERS
# damping ratio
damp = 0.02
# lower frequency
omega1 = 2 * np.pi * 0.2
# upper frequency
omega2 = 2 * np.pi * 20
# damping coefficients
a0 = 2*damp*omega1*omega2/(omega1 + omega2)
a1 = 2*damp/(omega1 + omega2)
print(f"FATTORI DI SMORZAMENTO: a_0 = {a0};  a_1 = {a1}")

#---DETERMINE STABLE ANALYSIS TIME STEP USING CFL CONDITION
# maximum shear wave velocity (m/s)
vsMax = 250.0
# duration of ground motion (s)
duration = motionDT*motionSteps
# minimum element size
minSize = 2

# trial analysis time step
kTrial = minSize / (vsMax ** 0.5)
# define time step and number of steps for analysis
if motionDT <= kTrial:
    nSteps = motionSteps
    dT = motionDT
else:
    nSteps = int(mm.floor(duration / kTrial) + 1)
    dT = duration / nSteps

set dT {dT}\n')
dT = {dT}\n')



gamma = 1.5
gamma1=0.5

beta  = (1/4)*(gamma+0.5)**2
beta1 = 0.25

ops.model("basicBuilder","-ndm",3,"-ndf",4)

ops.updateMaterialStage('-material', 1, '-stage', 0)
ops.updateMaterialStage('-material', 2, '-stage', 0)

##################################################### COSTRAINTS #################################################
ops.constraints('Penalty', 1.e18, 1.e18)

##################################################### TEST #######################################################
ops.test('NormDispIncr', 1.0e-6, 500, 1)

##################################################### ALGORITHMS #################################################
ops.algorithm('KrylovNewton')

##################################################### NUMBERER ###################################################
ops.numberer('Plain')

##################################################### SYSTEM #####################################################

ops.system('ProfileSPD')

##################################################### INTEGRATOR #################################################
ops.integrator('Newmark', gamma1, beta1)

##################################################### ANALYSIS ###################################################
ops.analysis('Transient')

##################################################### startT ###################################################
startT = tt.time()

##################################################### ANALYZE ###################################################
ops.analyze(1,1)

ops.updateMaterialStage('-material', 1, '-stage', 1)
ops.updateMaterialStage('-material', 2, '-stage', 1)

updateMaterialStage -material 2 -stage 1\n')

eleTags = []
getParamTags=[]


for i in ops.getParamTags():
    getParamTags.append(i)

LT = len(getParamTags)
if getParamTags == []:
    tag = 1
else:
    tag = int(getParamTags[LT-1])+1

for i in ops.getEleTags():
    eleTags.append(int(i))

LT2=len(eleTags)
for j in range(0,LT2-1): #ultimo elemento escluso perchè zeroLenght
    i = eleTags[j]
    ops.parameter(int(tag), 'element', i,'hPerm')
    ops.parameter(int(tag+1), 'element', i,'vPerm')

    ops.updateParameter(int(tag),0.001)
    ops.updateParameter(int(tag+1),0.001)
    tag = tag + 2

ops.timeSeries('Constant',1)
ops.pattern('Plain',1,1,'-fact',1.0)
ops.reactions()

g2o.get_physical_groups_map(gmsh.model) #dictionary - numero tag e dimensione ph ent.
TagPh = PhysGr['FixX'][1]
DimPh = PhysGr['FixX'][0]
entities_x=gmsh.model.getEntitiesForPhysicalGroup(DimPh,TagPh)
surfTags = gmsh.model.mesh.getElements(DimPh, entities_x[0]) # punti della ph ent selezionata
entTagDown,elTagsDown,nodeTagsDown=gmsh.model.mesh.getElements(2, entities_x[0])

xminSurf,yminSurf,zminSurf,xmaxSurf,ymaxSurf,zmaxSurf=gmsh.model.get_bounding_box(2,entities_x[0])

surfNodeTags=[]
surfNodeCorn = []
surfNodeEdge = []
NodeListX = []
mem=[]
n = 0
i = 0
j = 0

len1 = len(surfTags[2][0])
for i in range(0, len1):
        n = surfTags[2][0][i]
        surfNodeTags.append(n)
        NodeListX.append(n)
        mem.append(n)
        j=j+1
        if j == 4:
            surfNodeCorn.append(mem)
            mem=[]
        elif j == 8:
            surfNodeEdge.append(mem)
            j = 0
            mem=[]

LED=len(surfNodeEdge)
vecEdge=[]
for j in range(0,LED):
    for i in surfNodeEdge[j]:
        vecEdge.append(int(i))

vecEdgeA = np.array(vecEdge, dtype=int).reshape(-1)
vecEdgeU = np.unique(vecEdgeA)

reacV1={}
mem = 0
for i in vecEdgeU:

    mem = ops.nodeReaction(int(i), -1)  
    reacV1[int(i)] = mem
    ops.load(int(i), mem[0],0.0,0.0)  

    ops.remove('sp',int(i),1)




print("superficie down vettori creati")
xminSurf,yminSurf,zminSurf,xmaxSurf,ymaxSurf,zmaxSurf=gmsh.model.get_bounding_box(2,entities_down[0])
print('box nodi di base')
print(xminSurf,yminSurf,zminSurf,xmaxSurf,ymaxSurf,zmaxSurf)


surfNodeCornTags=[]
LLA=len(surfNodeCorn)
LLB=len(surfNodeCorn[0])
for i in range(0,LLA):
    for j in range(0,LLB):
        surfNodeCornTags.append(surfNodeCorn[i][j])

surfNodeCornA = np.array(surfNodeCornTags,dtype=int).reshape(-1)
surfNodeCornU = np.unique(surfNodeCornA)
xCoordCorn = 0
yCoordCorn = 0
zCoordCorn = 0
Coord_Corn_l = []
CoordCorn={}
len2=0
nodi_bordo1=[]
nodi_bordo2=[]
nodi_bordo3=[]
nodi_bordo4=[]
i = 0
j = 0

reacV2={}
mem=0
for j in surfNodeCornU:

        xCoordCorn = gmsh.model.mesh.getNode(j)[0][0]
        yCoordCorn = gmsh.model.mesh.getNode(j)[0][1]
        zCoordCorn = gmsh.model.mesh.getNode(j)[0][2]
        Coord_Corn_l = [xCoordCorn, yCoordCorn, zCoordCorn]
        CoordCorn[j] = Coord_Corn_l
        if (xCoordCorn == xmaxSurf) or (xCoordCorn == xmaxSurf and (yCoordCorn == ymaxSurf or yCoordCorn == yminSurf )):

           mem = ops.nodeReaction(int(j), -1)  
           ops.load(int(j), mem[0],0.0,0.0,0.0)  
           reacV2[int(j)] = mem

           ops.remove('sp', int(j), 1)

           

        elif (xCoordCorn == xminSurf) or (xCoordCorn == xminSurf and (yCoordCorn == ymaxSurf or yCoordCorn == yminSurf )):

           mem = ops.nodeReaction(int(j), -1)  
           ops.load(int(j), mem[0],0.0,0.0,0.0)  
           reacV2[int(j)] = mem

           ops.remove('sp', int(j), 1)

xCoord = 0
yCoord = 0
zCoord = 0
Coord_l=[]
Coord={}
reacV3={}
for i in vecEdgeU:
    xCoord = gmsh.model.mesh.getNode(i)[0][0]
    yCoord = gmsh.model.mesh.getNode(i)[0][1]
    zCoord = gmsh.model.mesh.getNode(i)[0][2]
    Coord_l = [xCoord, yCoord, zCoord]
    Coord[i] = Coord_l

    mem = ops.nodeReaction(int(i), -1)  
    ops.load(int(i), mem[0],0.0,0.0)  
    reacV3[int(i)] = mem

    ops.remove('sp',int(i), 1)

nodi_SurfxMax = []
nodi_SurfxMin = []
for i in NodeListX:
    xCoord = gmsh.model.mesh.getNode(i)[0][0]
    if xCoord == xmaxSurf:
        nodi_SurfxMax.append(i)
    elif xCoord == gmsh.model.mesh.getNode(i)[0][0]:
        nodi_SurfxMin.append(i)

i = 0
j = 0
k = 0
DOF = [1, 2]
eqDofDic = {}
for i in nodi_SurfxMax:
    for j in nodi_SurfxMin:
        yCoord_xmaxSurf = gmsh.model.mesh.getNode(i)[0][1]
        yCoord_xminSurf = gmsh.model.mesh.getNode(j)[0][1]
        zCoord_xmaxSurf = gmsh.model.mesh.getNode(i)[0][2]
        zCoord_xminSurf = gmsh.model.mesh.getNode(j)[0][2]
        if (yCoord_xmaxSurf == yCoord_xminSurf) and (zCoord_xmaxSurf == zCoord_xminSurf):
          
            eqDofDic[i] = j
            ops.equalDOF(int(i), int(j), *DOF)

cFactor = colArea * dashpotCoeff


# reset time and analysis
ops.setTime(0.0)
ops.wipeAnalysis()
ops.remove('recorders')


# recorder time step
recDT = 10*motionDT

tcl = open('cubotto.tcl', 'a')
tcl.write('set recDT [expr 10*$motionDT]\n')
tcl.close

# record nodal displacment, acceleration, and porepressure
ops.recorder('Node','-file','displacement.out','-time', '-dT',recDT,'-node',*nodeList,'-dof', 1, 2,3, 'disp')
ops.recorder('Node','-file','acceleration.out','-time', '-dT',recDT,'-node',*nodeList,'-dof', 1, 2,3, 'accel')
ops.recorder('Node','-file','porePressure.out','-time', '-dT',recDT,'-node',*nodeList,'-dof', 4, 'vel')



# record elemental stress and strain (files are names to reflect GiD gp numbering)
ops.recorder('Element','-file','stress.txt','-time', '-dT',recDT,'-eleRange', 211,335,'material','1','stress')
ops.recorder('Element','-file','strain.txt','-time', '-dT',recDT,'-eleRange', 211,335,'material','1','strain')
ops.recorder('Element','-file','strain.txt','-time', '-dT',recDT,'-eleRange', 211,335,'material','1','plastic')

ops.recorder('Element','-file','stress.out','-time', '-dT',recDT,'-eleRange', 211,335,'material','1','stress')
ops.recorder('Element','-file','strain.out','-time', '-dT',recDT,'-eleRange', 211,335,'material','1','strain')
ops.recorder('Element','-file','plastic.out','-time', '-dT',recDT,'-eleRange', 211,335,'material','1','plastic')


ops.model('basic', '-ndm', 3, '-ndf', 4)


# define constant scaling factor for applied velocity
cFactor = colArea * dashpotCoeff

# define velocity time history file
velocityFile = 'yerbaNSvelocity';
data_gm = np.loadtxt('yerbaNSvelocity.out')
# motionSteps=len(data_gm)
# print('Number of point for GM:',motionSteps)
velocityfileW = 'yerbaNSvelocity.out'
# timeseries object for force history
ops.timeSeries('Path', 2, '-dt', motionDT, '-filePath', velocityFile +'.out', '-factor', cFactor)
ops.pattern('Plain', 10, 2)
ops.load(458, 1.0,0.0, 0.0)

print("Dynamic loading created...")

ops.constraints('Penalty', 1.0E16, 1.0E16)
ops.test('NormDispIncr', 1e-3, 100, 1)
ops.algorithm('KrylovNewton')
ops.numberer('RCM')
ops.system('ProfileSPD')
ops.integrator('Newmark', gamma, beta)
ops.rayleigh(a0, a1, 0.0, 0.0)
ops.analysis('Transient')



# perform analysis with timestep reduction loop
ok = ops.analyze(nSteps, dT)

# if analysis fails, reduce timestep and continue with analysis
if ok != 0:
    print("did not converge, reducing time step")
    curTime = ops.getTime()
    mTime = curTime
    print("curTime: ", curTime)
    curStep = curTime / dT
    print("curStep: ", curStep)
    rStep = (nSteps - curStep) * 2.0
    remStep = int((nSteps - curStep) * 2.0)
    print("remStep: ", remStep)
    dT = dT / 2.0
    print("dT: ", dT)

    ok = ops.analyze(remStep, dT)
    # if analysis fails again, reduce timestep and continue with analysis
    if ok != 0:
        print("did not converge, reducing time step")
        curTime = ops.getTime()
        print("curTime: ", curTime)
        curStep = (curTime - mTime) / dT
        print("curStep: ", curStep)
        remStep = int((rStep - curStep) * 2.0)
        print("remStep: ", remStep)
        dT = dT / 2.0
        print("dT: ", dT)

        ok = ops.analyze(remStep, dT)

endT = tt.time()
print("Finished with dynamic analysis...")
print("Analysis execution time: ", (endT - startT))
ops.wipe()


