#!/usr/bin/env python

#> \file 
#> \author David Ladd, Reused: Hashem Yousefi 
#> \brief This is an example to use linear fitting to fit the beginning of linear heart tube.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>
### to be used for segmenting embryonic heart and fitting with an initial meshes.
#<

import sys, os
import exfile
import numpy
from numpy import linalg
import math
import random

# Intialise OpenCMISS/iron 
from opencmiss.iron import iron

# defining the output file to be written in the ExDataFile
def writeExdataFile(filename,dataPointLocations,offset):
    "Writes data points to an exdata file"

    numberOfDimensions = dataPointLocations[1].shape[0]
    try:
        f = open(filename,"w")    
        header = '''Group name: DataPoints
 #Fields=1
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  1.  Value index=1, #Derivatives=0, #Versions=1
  2.  Value index=2, #Derivatives=0, #Versions=1
'''
        if numberOfDimensions == 3:
            header+= '''  3.  Value index=3, #Derivatives=0, #Versions=1
'''
        f.write(header)

        numberOfDataPoints = len(dataPointLocations)
        for i in range(numberOfDataPoints):
            line = " Node: " + str(offset+i+1) + '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataPointLocations[i,j]) + '\n'
                f.write(line)
        f.close()
            
    except IOError:
        print ('Could not open file: ' + filename)

#=================================================================
# Control Panel
#=================================================================
# set the number of elements and the number of nodes for the cylinder 
numberOfDimensions = 3
numberOfGaussXi = 3 
numberOfCircumfrentialElementsPerQuarter = 1
numberOfCircumfrentialElements = 4*numberOfCircumfrentialElementsPerQuarter
numberOfCircumfrentialNodes = numberOfCircumfrentialElements
numberOfLengthElements = 5
numberOfLengthNodes = numberOfLengthElements+1
numberOfWallElements = 1
numberOfWallNodes = numberOfWallElements+1
origin = [0.0,0.0,0.0]
meshOrigin = [0.0,0.0,0.0]
print "mesh resolution and parameters fixed"

# The number of data points which are digitised from the heart segments 
# fix interior nodes so that fitting only applies to surface
# If start iteration > 1, read in geometry from a previous fit iteration
numberOfDataPoints = 1024
numberOfIterations = 1
fixInterior = True
zeroTolerance = 0.00001
hermite = True
# Set Sobolev smoothing parameters
tau = 0.5
kappa = 0.1
iteration = 1
if iteration > 1:
    exfileMesh = True
    exnode = exfile.Exnode("DeformedGeometry" + str(iteration-1) + ".part0.exnode")
    exelem = exfile.Exelem("UndeformedGeometry.part0.exelem")
else:
    exfileMesh = False
print "other parameters setted up "

#=================================================================
# 
#=================================================================
(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    independentFieldUserNumber,
    dataPointFieldUserNumber,
    materialFieldUserNumber,
    analyticFieldUserNumber,
    dependentDataFieldUserNumber,
    dataProjectionUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,18)

# Get the computational nodes information
print dir(iron),'\n\n'
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a RC CS
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "FittingRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# define a basis 
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = numberOfDimensions
if hermite:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
else:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
basis.CreateFinish()
print "CS, Region and basis setted up"

#=================================================================
# Mesh
#=================================================================
# creating the number of elements and the mesh origins ... and/or
# Start the creation of a manually generated mesh in the region
numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements

print "numberOfElements = ", numberOfElements
print "numberOfNodes = ", numberOfNodes


if (exfileMesh):
    # Read previous mesh
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, numberOfDimensions)
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(exelem.num_elements)
    # Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region, exnode.num_nodes)
    nodes.CreateFinish()
    # Define elements for the mesh
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    for elem in exelem.elements:
        elements.NodesSet(elem.number, elem.nodes)
    elements.CreateFinish()
    mesh.CreateFinish()
else:
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber,region,3)
    mesh.origin = meshOrigin
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(numberOfElements)
# Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region,numberOfNodes)
    nodes.CreateFinish()
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    elementNumber = 0
    for wallElementIdx in range(1,numberOfWallElements+1):
       for lengthElementIdx in range(1,numberOfLengthElements+1):
            for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
                elementNumber = elementNumber + 1
                localNode1 = circumfrentialElementIdx + (lengthElementIdx - 1)*numberOfCircumfrentialElements + \
                    (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                if circumfrentialElementIdx == numberOfCircumfrentialElements:
                    localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                        (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                else: 
                    localNode2 = localNode1 + 1
                localNode3 = localNode1 + numberOfCircumfrentialNodes
                localNode4 = localNode2 + numberOfCircumfrentialNodes
                localNode5 = localNode1 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode6 = localNode2 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode7 = localNode3 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode8 = localNode4 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
#		print "Element Number = ",elementNumber
#         	print "Node numbers of the element", localNode1, localNode2, localNode3, localNode4, localNode5, localNode6, localNode7, localNode8 
                elements.NodesSet(elementNumber,localNodes)  
    elements.CreateFinish()
    mesh.CreateFinish() 


# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

print "mesh decomposition finished"

#=================================================================
# Geometric Field
#=================================================================
# the location of  nodes for the mesh  
manualNodePoints = numpy.zeros((6,4,3,2))

manualNodePoints[0,0,:,0] = [528,315,-5]
manualNodePoints[0,1,:,0] = [513,322,-5]
manualNodePoints[0,2,:,0] = [497,315,-5]
manualNodePoints[0,3,:,0] = [512,305,-5]

manualNodePoints[1,0,:,0] = [548,300,115]
manualNodePoints[1,1,:,0] = [526,324,115]
manualNodePoints[1,2,:,0] = [504,301,115]
manualNodePoints[1,3,:,0] = [526,287,115]

manualNodePoints[2,0,:,0] = [595,260,235]
manualNodePoints[2,1,:,0] = [550,290,235]
manualNodePoints[2,2,:,0] = [467,265,235]
manualNodePoints[2,3,:,0] = [550,250,235]

manualNodePoints[3,0,:,0] = [625,225,365]
manualNodePoints[3,1,:,0] = [550,225,365]
manualNodePoints[3,2,:,0] = [425,250,365]
manualNodePoints[3,3,:,0] = [550,188,365]

manualNodePoints[4,0,:,0] = [675,260,515]
manualNodePoints[4,1,:,0] = [590,260,515]
manualNodePoints[4,2,:,0] = [505,260,515]
manualNodePoints[4,3,:,0] = [590,220,515]

manualNodePoints[5,0,:,0] = [607,304,635]
manualNodePoints[5,1,:,0] = [587,332,635]
manualNodePoints[5,2,:,0] = [550,304,635]
manualNodePoints[5,3,:,0] = [587,288,635]

# node locations of the outer surface ... 
manualNodePoints[0,0,:,1] = [595,311,-5]
manualNodePoints[0,1,:,1] = [515,338,-5]
manualNodePoints[0,2,:,1] = [440,311,-5]
manualNodePoints[0,3,:,1] = [515,280,-5]

manualNodePoints[1,0,:,1] = [660,280,105]
manualNodePoints[1,1,:,1] = [515,375,105]
manualNodePoints[1,2,:,1] = [330,280,105]
manualNodePoints[1,3,:,1] = [545,220,105]

manualNodePoints[2,0,:,1] = [710,250,225]
manualNodePoints[2,1,:,1] = [550,370,225]
manualNodePoints[2,2,:,1] = [325,275,225]
manualNodePoints[2,3,:,1] = [550,178,225]

manualNodePoints[3,0,:,1] = [763,235,375]
manualNodePoints[3,1,:,1] = [550,310,375]
manualNodePoints[3,2,:,1] = [325,257,375]
manualNodePoints[3,3,:,1] = [552,140,375]

manualNodePoints[4,0,:,1] = [830,250,525]
manualNodePoints[4,1,:,1] = [510,358,525]
manualNodePoints[4,2,:,1] = [298,257,525]
manualNodePoints[4,3,:,1] = [540,145,525]

manualNodePoints[5,0,:,1] = [800,300,645]
manualNodePoints[5,1,:,1] = [575,375,645]
manualNodePoints[5,2,:,1] = [360,270,645]
manualNodePoints[5,3,:,1] = [575,180,645]

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
for dimension in range(3):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,dimension+1,1)
geometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
geometricField.CreateFinish()

# Create the geometric field
for wallNodeIdx in range(1,numberOfWallNodes+1):
    for lengthNodeIdx in range(1,numberOfLengthNodes+1):
        for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
            nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + \
                (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            theta = float(circumfrentialNodeIdx-1)/float(numberOfCircumfrentialNodes)*2.0*math.pi
            if (nodeNumber < 25):
                x = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, 0]
                y = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, 0]
                z = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, 0]
            else: 
                x = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, 1]
                y = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, 1]
                z = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, 1]
            xtangent = -math.sin(theta) 
            ytangent = math.cos(theta)
            xnormal = math.cos(theta)
            ynormal = math.sin(theta)        
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,1,nodeNumber,1,x)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,1,nodeNumber,2,y)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,1,nodeNumber,3,z)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,1.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,0.0)


# Update the geometric field
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Get nodes
nodes = iron.Nodes()
region.NodesGet(nodes)
numberOfNodes = nodes.numberOfNodes

# Get or calculate geometric parameters
if (exfileMesh):
    # Read the geometric field from the exnode file
    geometricField.ParameterSetUpdateStart(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    for node_num in range(1, exnode.num_nodes + 1):
        version = 1
        derivative = 1
        for component in range(1, numberOfDimensions + 1):
            component_name = ["x", "y", "z"][component - 1]
            value = exnode.node_value("Coordinate", component_name, node_num, derivative)
            geometricField.ParameterSetUpdateNode(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    version, derivative, node_num, component, value)
    geometricField.ParameterSetUpdateFinish(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
else:
    # Create undeformed geometry from the generated mesh
#    mesh.GeometricParametersCalculate(geometricField)
    # Export undeformed mesh geometry
    print("Writing undeformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("UndeformedGeometry","FORTRAN")
    fields.ElementsExport("UndeformedGeometry","FORTRAN")
    fields.Finalise()

#=================================================================
# Data Points
#=================================================================
# Create the data points
dataPoints = iron.DataPoints()
dataPoints.CreateStart(region,numberOfDataPoints)
dataPointLocations = numpy.zeros((numberOfDataPoints,3))
print("Number of data points: " + str(numberOfDataPoints))
# reading from a text file containing the point clouds  
with open("FullPoints8dpc.txt", "r") as ins:
	arrayOfInputData = []
	for line in ins:
		arrayOfInputData.append(line)
x = 0
y = 0
z = 0
for i in range (numberOfDataPoints):
	for j in range (5):
		sample = arrayOfInputData[i*5 + j]
		if (math.fmod(j,5) == 1):
			x = float (sample[12:25])				
		elif (math.fmod(j,5) == 2):
			y = float (sample[12:25])
		elif (math.fmod(j,5) == 3):
			z = float (sample[12:17])
		dataPointLocations[i,:] = [x,y,z]
# Set up data points with geometric values
for dataPoint in range(numberOfDataPoints):
    dataPointId = dataPoint + 1
    dataList = dataPointLocations[dataPoint,:]
    dataPoints.ValuesSet(dataPointId,dataList)
dataPoints.CreateFinish()
 
# write data points to exdata file for CMGUI
offset = 0
writeExdataFile("DataPoints.part"+str(computationalNodeNumber)+".exdata",dataPointLocations,offset)

#=================================================================
# Data Projection on Geometric Field
#=================================================================
print("Projecting data points onto geometric field")
# Set up data projection
dataProjection = iron.DataProjection()
dataProjection.CreateStart(dataProjectionUserNumber,dataPoints,mesh)
dataProjection.projectionType = iron.DataProjectionProjectionTypes.ALL_ELEMENTS
dataProjection.CreateFinish()

# Evaluate data projection based on geometric field
dataProjection.DataPointsProjectionEvaluate(geometricField)
# Create mesh topology for data projection
mesh.TopologyDataPointsCalculateProjection(dataProjection)
# Create decomposition topology for data projection
decomposition.TopologyDataProjectionCalculate()
print("Projection complete")
 
#=================================================================
# Equations Set
#=================================================================
# Create vector fitting equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                             iron.EquationsSetTypes.DATA_FITTING_EQUATION,
                             iron.EquationsSetSubtypes.DATA_POINT_VECTOR_STATIC_FITTING]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

#=================================================================
# Dependent Field
#=================================================================
# Create dependent field (will be deformed fitted values based on data point locations)
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
equationsSet.DependentCreateFinish()
# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)

# Initialise dependent field to undeformed geometric field
for component in range (1,numberOfDimensions+1):
    geometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            component, dependentField, iron.FieldVariableTypes.U,
                                                            iron.FieldParameterSetTypes.VALUES, component)

#=================================================================
# Independent Field
#=================================================================
# Create data point field (independent field, with vector values stored at the data points)
independentField = iron.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber,independentField)
independentField.VariableLabelSet(iron.FieldVariableTypes.U,"data point vector")
independentField.VariableLabelSet(iron.FieldVariableTypes.V,"data point weight")
independentField.DataProjectionSet(dataProjection)
equationsSet.IndependentCreateFinish()
# Initialise data point vector field to 0
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
# Initialise data point weight field to 1
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,1.0)
# loop over each element's data points and set independent field values to data point locations on surface of the sphere
for element in range(numberOfElements):
    elementId = element + 1
    elementDomain = decomposition.ElementDomainGet(elementId)
    if (elementDomain == computationalNodeNumber):
        numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementId)
        for dataPoint in range(numberOfProjectedDataPoints):
            dataPointId = dataPoint + 1
            dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementId,dataPointId)
            dataList = dataPoints.ValuesGet(dataPointNumber,3)
            # set data point field values
            for component in range(numberOfDimensions):
                componentId = component + 1
                dataPointNumberIndex = dataPointNumber - 1
                value = dataList[component]
                independentField.ParameterSetUpdateElementDataPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementId,dataPointId,componentId,value)

#=================================================================
# Material Field
#=================================================================
# Create material field (Sobolev parameters)
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Smoothing Parameters")
equationsSet.MaterialsCreateFinish()
# Set kappa and tau - Sobolev smoothing parameters
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)

#=================================================================
# Equations
#=================================================================
# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.FULL
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#=================================================================
# Problem setup
#=================================================================
# Create fitting problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FITTING,
                        iron.ProblemTypes.DATA_FITTING,
                        iron.ProblemSubtypes.DATA_POINT_VECTOR_STATIC_FITTING]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.NONE # NONE / MATRIX
#solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearType = iron.LinearSolverTypes.DIRECT
#solver.LibraryTypeSet(iron.SolverLibraries.UMFPACK) # UMFPACK/SUPERLU
solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
#solver.linearIterativeAbsoluteTolerance = 1.0E-10
#solver.linearIterativeRelativeTolerance = 1.0E-05
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
#solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.FULL
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#=================================================================
# Boundary Conditions
#=================================================================
# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

version = 1
meshComponent = decomposition.MeshComponentGet()
# Fix the interior nodes- use to only apply fit to surface nodes
if (fixInterior):
    # first find which nodes are non-surface nodes
    for node in range(numberOfNodes):
        nodeId = node + 1
        nodeDomain = decomposition.NodeDomainGet(nodeId,meshComponent)
        if (nodeDomain == computationalNodeNumber):
            geometricValue = numpy.zeros((numberOfDimensions))
            for component in range(numberOfDimensions):
                componentId = component+1
                geometricValue[component]=geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                               iron.FieldParameterSetTypes.VALUES,
                                                                               1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                                               nodeId,componentId)

            derivList = []
            deriv = 0
            if hermite:
                if numberOfDimensions ==3:
                    if (abs(geometricValue[0]) < (abs(meshOrigin[0]) - zeroTolerance)) and (abs(geometricValue[1]) < (abs(meshOrigin[1]) - zeroTolerance)) and (abs(geometricValue[2]) < (abs(meshOrigin[2]) - zeroTolerance)):
                        # Interior nodes
                        derivList = [1,2,3,4,5,6,7,8]
                    # Radial nodes
                    elif abs(geometricValue[0]) < zeroTolerance and abs(geometricValue[1]) < zeroTolerance:
                        deriv = 5
                    elif abs(geometricValue[1]) < zeroTolerance and abs(geometricValue[2]) < zeroTolerance:
                        deriv = 2
                    elif abs(geometricValue[0]) < zeroTolerance and abs(geometricValue[2]) < zeroTolerance:
                        deriv = 3

                    if deriv > 0 and deriv not in derivList:
                        derivList.append(deriv)

                elif numberOfDimensions ==2:
                    if (abs(geometricValue[0]) < (abs(meshOrigin[0]) - zeroTolerance)) and (abs(geometricValue[1]) < (abs(meshOrigin[1]) - zeroTolerance)):
                        # Interior nodes
                        derivList = [1,2,3,4]
                    # Radial nodes
                    elif abs(geometricValue[1]) < zeroTolerance:
                        deriv = 2
                    elif abs(geometricValue[0]) < zeroTolerance:
                        deriv = 3

                    if deriv > 0 and deriv not in derivList:
                        derivList.append(deriv)

            elif numberOfDimensions == 3:
                if (abs(geometricValue[0]) < (abs(meshOrigin[0]) - zeroTolerance)) and (abs(geometricValue[1]) < (abs(meshOrigin[1]) - zeroTolerance)) and (abs(geometricValue[2]) < (abs(meshOrigin[2]) - zeroTolerance)):
                    derivList = [1]
            elif numberOfDimensions == 2:
                if (abs(geometricValue[0]) < (abs(meshOrigin[0]) - zeroTolerance)) and (abs(geometricValue[1]) < (abs(meshOrigin[1]) - zeroTolerance)):
                    derivList = [1]

            for globalDeriv in derivList: 
                for component in range(1,numberOfDimensions+1):
                    value=geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                               iron.FieldParameterSetTypes.VALUES,
                                                               version,globalDeriv,nodeId,component)
                    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,
                                               version,globalDeriv,nodeId,component,
                                               iron.BoundaryConditionsTypes.FIXED,value)

solverEquations.BoundaryConditionsCreateFinish()

#=================================================================
# S o l v e    a n d    E x p o r t    D a t a
#=================================================================
for iteration in range (1,numberOfIterations+1):
    # Solve the problem
    print("Solving fitting problem, iteration: " + str(iteration))
    problem.Solve()
    # Copy dependent field to geometric 
    for component in range(1,numberOfDimensions+1):
        dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,component,geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,component)
    # Export fields
    print("Writing deformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.ElementsExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.Finalise()
#-----------------------------------------------------------------
iron.Finalise()
