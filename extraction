# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import*
from abaqusConstants import*
import __main__
from os import listdir
from os import system
from odbAccess import openOdb
defOdb_files(path):
OdbFilesList =[]
graphcount =1
for name in listdir(path):
if name.endswith('.odb'):
OdbFilesList.append(name)
OdbFilesList.sort()
for current_file in OdbFilesList:
directory(path,graphcount,current_file)
graphcount +=1
defdirectory(path,graphcount,fname):
address =""
address = path + fname
o1 = session.openOdb(name=address)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
Extraction(address,graphcount,fname)
odb = session.odbs[address]
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.odbs[address].close()
defExtraction(address,graphcount,fname):
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
odb = session.odbs[address]
session.xyDataListFromField(odb=odb, outputPosition=NODAL,
variable=(('RF',
NODAL,((COMPONENT,'RF2'),)),), nodeSets=('SUPPORTBC',)) odb = session.odbs[address]
session.xyDataListFromField(odb=odb, outputPosition=NODAL,
variable=(('U',
NODAL,((COMPONENT,'U2'),)),), nodeSets=(
'REFERENCE_POINT_ 1',))
Total_force =0
Node_counter =0
Values_List =[]
for key, value in session.xyDataObjects.items():
Total_force += value
Node_counter +=1
MaximumPower = max(Total_force, key=lambda x: x[1])
lst = list(MaximumPower)
MaximumValue = lst[1]/1000
MaximumValue = str(MaximumValue)
print(fname +" = "+ MaximumValue)
for delkey, delvalue in session.xyDataObjects.items():
if delkey in session.xyDataObjects:
del session.xyDataObjects[delkey]
del MaximumPower
del MaximumValue
path='G:/Amir/Table 1/'
Odb_files(path)
