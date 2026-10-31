# -*- coding: utf-8 -*-

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
from array import *
import math
import numpy
import os        # Operating system
import shutil    # copying or moving files
from abaqusConstants import *



# ---------------------- DEFINE PARAMETERS ------------------------ #

s_f = 0.0001 # scale factor (m --> 0.1mm)

# INCLUSION RADIUS
global RAD;
RAD = 10.0 * s_f

# INCLUSION/EDGE SEPARATION
EDGE_CENTRE_SEP = 12.0 * s_f # min inclusion distance from edge (should be bigger than RADius)

# INCLUSION SEPARATION
INCLUSION_SEP = 2.2  # number of RADii which inclusion centres are sep by (i.e. > 2)

# MATRIX DIMENSIONS
X_LENGTH = 600 * s_f
Y_LENGTH = 480 * s_f
Z_LENGTH = 65 * s_f

# LOCATE MATRIX
MATRIX_POS = 0.1 * s_f # point at which will definitely select matrix

# LOCATE SPHERE
SPHERE_CELLS_POS = 0.2 * s_f # point at which will select sphere cells from sphere origin (i.e. < RAD)

# NUMBER OF MODELS
MAX_ITERATIONS = 1 

# NUMBER OF INCLUSIONS
MAX_INCL = 10

# THERMAL PROPERTIES:
# i) MATRIX
MATRIX_CONDUC = 0.2 # W/m/K
MATRIX_SP_HEAT = 1200.0  # J/kg/K
MATRIX_DENSITY = 1200.0 # kg/m^3
# ii) INCLUSIONS
INC_CONDUC = 2000.0 # W/m/K
INC_SP_HEAT = 700.0  # J/kg/K
INC_DENSITY = 2200.0  # kg/m^3

# PULSE PROPERTIES:
MAGNITUDE = 2500.0 # W/m^2
PULSE_DURATION = 2.0 # seconds
INCREMENT = 0.5 # seconds

# SIMULATION DURATION
TOTAL_DURATION = 240.0 # seconds 	

# SEED SIZE
SEED_SIZE = X_LENGTH/80

# CONVECTION BCs
SURR_TEMP = 20.0 # temp of surroundings in degrees
FILM_COEFF = 10.0 # film coefficient for convection

# ------------------------ PARTITION FUNCTION ----------------------------- #

def partition(i,q):
    """ Create a spherical partition """

    # defining a and b (not sure what for yet)
    a=RAD*sin(pi/6)
    b=RAD*cos(pi/6)

    # Creating Datum Planes: used to draw spheres away from origin
    # line at which the planes intersect is used to revolve semi-circle into sphere
    dp1=mdb.models['Model-%d' %(q)].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=x_coordinate[i], principalPlane=YZPLANE)
    dp2=mdb.models['Model-%d' %(q)].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=y_coordinate[i], principalPlane=XZPLANE)
	
    # opens a sketch - uses transform (because its a sketch on a sketch)
    # sketch plane is datum plane d1
    # define "up" using edge of cube at 0,0,Z_LENGTH/4
    # sets the origin of the datum plane at where we want the centre of the sphere
    mdb.models['Model-%d' %(q)].ConstrainedSketch(gridSpacing=2.0*s_f, name='__profile__', sheetSize=25.0*s_f, transform=
        mdb.models['Model-%d' %(q)].parts['Part-1'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].datums[dp1.id], sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-%d' %(q)].parts['Part-1'].edges.findAt((0.0, 0.0, Z_LENGTH/4), ), sketchOrientation=TOP, origin=(x_coordinate[i],y_coordinate[i],z_coordinate[i])))

    # projects existing model edges onto sketch
    mdb.models['Model-%d' %(q)].parts['Part-1'].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])

    # draws quarter circle arc about datum plane origin
    c_1=mdb.models['Model-%d' %(q)].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(0,RAD), point2=(RAD,0))

    # draws another connected quarter circle arc about datum plane origin
    c_2=mdb.models['Model-%d' %(q)].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(RAD,0), point2=(0,-RAD))

    # draws a line joining the arc ends
    mdb.models['Model-%d' %(q)].sketches['__profile__'].Line(point1=(0,RAD), point2=(0,-RAD))

    # partitions the matrix into a "half cyliner" (i.e. projects the semi circle perpendicular to plane its drawn in)
    # doesnt actually cut it out yet
    mdb.models['Model-%d' %(q)].parts['Part-1'].PartitionCellBySketch(cells=
        mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((MATRIX_POS, MATRIX_POS, MATRIX_POS), )), sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'], 
        sketchOrientation=TOP, sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].datums[dp1.id], sketchUpEdge=mdb.models['Model-%d' %(q)].parts['Part-1'].edges.findAt((0.0, 0.0, Z_LENGTH/4), ))
    
    # opens another sketch, this time on datum 2
    # exactly the same as other sketch just different plane (perp to other)
    mdb.models['Model-%d' %(q)].ConstrainedSketch(gridSpacing=2.0*s_f, name='__profile__', sheetSize=25.0*s_f,
        transform=mdb.models['Model-%d' %(q)].parts['Part-1'].MakeSketchTransform(sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].datums[dp2.id],
        sketchPlaneSide=SIDE1,sketchUpEdge=mdb.models['Model-%d' %(q)].parts['Part-1'].edges.findAt((0.0, 0.0, Z_LENGTH/4), ), sketchOrientation=TOP, origin=(x_coordinate[i], y_coordinate[i], z_coordinate[i])))

    # project existing edges onto sketch
    mdb.models['Model-%d' %(q)].parts['Part-1'].projectReferencesOntoSketch(filter=COPLANAR_EDGES,sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])
    c_3=mdb.models['Model-%d' %(q)].sketches['__profile__'].CircleByCenterPerimeter(center=(0, 0), point1=(-RAD, 0))

    # do the same semi circle partition thing (semi circle cylinder thing)
    mdb.models['Model-%d' %(q)].parts['Part-1'].PartitionCellBySketch(cells=mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((SPHERE_CELLS_POS,SPHERE_CELLS_POS, 
        SPHERE_CELLS_POS), )), sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'],sketchOrientation=TOP, sketchPlane=
        mdb.models['Model-%d' %(q)].parts['Part-1'].datums[dp2.id], sketchUpEdge=mdb.models['Model-%d' %(q)].parts['Part-1'].edges.findAt((0.0, 0.0, Z_LENGTH/4), ))
		
    # Creating the spherical partition
    # trust this step
    m= mdb.models['Model-%d' %(q)].parts['Part-1']
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((x_coordinate[i]+RAD, y_coordinate[i], z_coordinate[i]),),cells=m.cells.findAt((SPHERE_CELLS_POS, SPHERE_CELLS_POS, SPHERE_CELLS_POS),),edges=(m.edges.findAt((x_coordinate[i], y_coordinate[i]-a,z_coordinate[i]+b), ),))
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((x_coordinate[i]+RAD, y_coordinate[i], z_coordinate[i]),),cells=m.cells.findAt((SPHERE_CELLS_POS, SPHERE_CELLS_POS, SPHERE_CELLS_POS),),edges=(m.edges.findAt((x_coordinate[i], y_coordinate[i]+a,z_coordinate[i]+b), ),))	
    

# ---------------- END OF PARTITION FUNCTION ------------------ #


dis=numpy.zeros(1000)

for q in range (1,MAX_ITERATIONS+1):

    # create a model
    # index with q
    mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-%d' %(q))    
    
    # create a matrix
    mdb.models['Model-%d' %(q)].ConstrainedSketch(name='__profile__', sheetSize=20.0)
    mdb.models['Model-%d' %(q)].sketches['__profile__'].sketchOptions.setValues(
        decimalPlaces=4)
    mdb.models['Model-%d' %(q)].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
        point2=(X_LENGTH, Y_LENGTH))
    mdb.models['Model-%d' %(q)].Part(dimensionality=THREE_D, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-%d' %(q)].parts['Part-1'].BaseSolidExtrude(depth=Z_LENGTH, sketch=
        mdb.models['Model-%d' %(q)].sketches['__profile__'])
    del mdb.models['Model-%d' %(q)].sketches['__profile__']        

    # set number of inclusions and empty arrays storing centres
    num_incl = 0
    x_coordinate = []
    y_coordinate = []
    z_coordinate = []

    # whilst still not the max num of inclusions
    while (num_incl < MAX_INCL):

        # create a random x y and z within the volume
        random_x=random.uniform(EDGE_CENTRE_SEP, X_LENGTH-EDGE_CENTRE_SEP)
        random_y=random.uniform(EDGE_CENTRE_SEP, Y_LENGTH-EDGE_CENTRE_SEP)
        random_z=random.uniform(EDGE_CENTRE_SEP, Z_LENGTH-EDGE_CENTRE_SEP)

	# set the checker for if the centre point is intersecting initially to false
        isPointIntersecting = False

	# calculate the distances to all of the centres (from the new random centre)
        for j in range (0,len(x_coordinate)):
    
	    # distance caluclation (pythag 3D)
            dis[j]=sqrt((random_x-x_coordinate[j])**2+(random_y-y_coordinate[j])**2+(random_z-z_coordinate[j])**2)

            # if the distance is less than "INCLUSION_SEP" of the RADii (i.e. almost intersecting)        
            if dis[j] < (INCLUSION_SEP*RAD):
		
		# set that its intersecting and break the loop
	        # breaking loop creates new point, tries again
                isPointIntersecting = True
                break
	
	# if the point is not intersecting
        if (isPointIntersecting == False):

            # store the centre point
            x_coordinate.append(random_x)
            y_coordinate.append(random_y)
            z_coordinate.append(random_z)

	    # count the inclusion
            num_incl = num_incl + 1        
    

    # for each inclusion
    for i in range (num_incl):

        # call the partition funciton to cut it out the matrix
        partition(i,q)

# ---------------------------- MATERIALS AND SECTION ASSIGNMENT ---------------------------- #

    # create a material named matrix
    mdb.models['Model-%d' %(q)].Material(name='matrix')

    # MATRIX THERMAL PROPERTIES
    # Assign thermal conductivity
    mdb.models['Model-%d' %(q)].materials['matrix'].Conductivity(table=((MATRIX_CONDUC,), ))

    # Assign specific heat
    mdb.models['Model-%d' %(q)].materials['matrix'].SpecificHeat(table=((MATRIX_SP_HEAT,), ))

    # Assign density
    mdb.models['Model-%d' %(q)].materials['matrix'].Density(table=((MATRIX_DENSITY,), ))
    
    # create a material for the inclusions
    mdb.models['Model-%d' %(q)].Material(name='GNP')

    # GNP thermal properties
    # Assign thermal conductivity
    mdb.models['Model-%d' %(q)].materials['GNP'].Conductivity(table=((INC_CONDUC,), ))

    # Assign specific heat
    mdb.models['Model-%d' %(q)].materials['GNP'].SpecificHeat(table=((INC_SP_HEAT,), ))

    # Assign density
    mdb.models['Model-%d' %(q)].materials['GNP'].Density(table=((INC_DENSITY,), ))
    
    # create a section for the matrix - material "matrix" 
    mdb.models['Model-%d' %(q)].HomogeneousSolidSection(material='matrix', name='Matrix', 
        thickness=None)
    # create a section for the inclusions - material "elastic"
    mdb.models['Model-%d' %(q)].HomogeneousSolidSection(material='GNP', name='Inclusion', 
        thickness=None)
        
    # assign the matrix section to the matrix
    mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        cells=mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((MATRIX_POS, MATRIX_POS, MATRIX_POS), ), )), sectionName='Matrix',
        thicknessAssignment=FROM_SECTION)
    
    # for each inclusion
    for i in range (num_incl):

	# assign the section to cells
        mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            cells=mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i], y_coordinate[i]-0.2*RAD, z_coordinate[i]+0.2*RAD), ), )),
            sectionName='Inclusion', thicknessAssignment=FROM_SECTION)

	# make sure we get all cells
        mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            cells=mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i], y_coordinate[i]+0.2*RAD, z_coordinate[i]+0.2*RAD), ), )),
            sectionName='Inclusion', thicknessAssignment=FROM_SECTION)

# ------------------------------ DEFINE ASSEMBLY AND STEP ---------------------------------- #

    # define shorthand
    block_model = mdb.models['Model-%d' %(q)]

    # creating an INSTANCE of the whole part
    # with the matrix and inclusions 
    mdb.models['Model-%d' %(q)].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-%d' %(q)].rootAssembly.Instance(dependent=ON, name='Composite_Instance', 
        part=mdb.models['Model-%d' %(q)].parts['Part-1'])

    # get important variable (without changing instance or assembly creation)
    block_assembly = block_model.rootAssembly
    block_instance = block_assembly.instances['Composite_Instance']

    # get all cells from the composite instance
    all_cells = block_instance.cells[:]

    # create a set that represents the whole block
    block_region = block_assembly.Set(name='whole_block', cells=all_cells)

    # Apply the initial temperature
    block_model.Temperature(
        name='InitialTemp',
        createStepName='Initial',
        region=block_region,
        distributionType=UNIFORM,
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
        magnitudes=(20.0,)
    )

    # Define the surface on the face at x=0 (front face of the block)
    target_face = block_instance.faces.findAt(((X_LENGTH/2, Y_LENGTH/2, Z_LENGTH),))
    block_assembly.Surface(name='heat_surface', side1Faces=target_face)

    # define a region from the surface to use later
    surface_region = block_assembly.surfaces['heat_surface']

    # Define an analysis step
    block_model.HeatTransferStep(
        name='heat_step',
        previous='Initial',
        timePeriod=TOTAL_DURATION,      # Simulate for number of seconds
        maxNumInc=5000,       # Max number of increments
        initialInc=INCREMENT,       # Initial time increment
        minInc=INCREMENT,
        maxInc=INCREMENT,
        deltmx=5.0  # Max temperature change per increment (tweak as needed)
    )
            
# ----------------------------- EDIT BOUNDARY CONDITIONS --------------------------------- #

    # Define a heat pulse that is active for only 0.1s
    block_model.TabularAmplitude(
        name='pulse_amplitude',
        timeSpan=STEP,  # amplitude applies over the step duration
        data=((0.0, 1.0), (PULSE_DURATION, 1.0), (PULSE_DURATION + PULSE_DURATION / 100, 0.0), (10.0, 0.0))
    )

    # Apply heat flux load (the heat pulse) to the defined surface in the heat_step
    block_model.SurfaceHeatFlux(
        name='heat_pulse_load',
        createStepName='heat_step',
        region=surface_region,
        magnitude=MAGNITUDE,          # peak flux in W/m^2, adjust as needed
        amplitude='pulse_amplitude'  # use the defined pulse amplitude
    )

    # Request field output
    block_model.FieldOutputRequest(
        name='F-Output-1',
        createStepName='heat_step',
        variables=('NT',),
        frequency=1  # Save temperature at every increment
    )

    # Define ambient temperature
    ambient_temp = SURR_TEMP  # 0°C in Kelvin
    film_coeff = FILM_COEFF      # W/mK (typical free air convection)

    # Get the front face at z = Z_LENGTH
    front_face = block_instance.faces.findAt(((X_LENGTH/2, Y_LENGTH/2, Z_LENGTH),))
    block_assembly.Surface(name='convection_front', side1Faces=front_face)
    region_front = block_assembly.surfaces['convection_front']

    # Apply convection boundary condition
    block_model.FilmCondition(
        name='ConvectionFront',
        createStepName='heat_step',
        surface=region_front,
        filmCoeff=film_coeff,
        sinkTemperature=ambient_temp,
        definition=EMBEDDED_COEFF,
    )

    # Back face at z = 0
    back_face = block_instance.faces.findAt(((X_LENGTH/2, Y_LENGTH/2, 0.0),))
    block_assembly.Surface(name='convection_back', side1Faces=back_face)
    region_back = block_assembly.surfaces['convection_back']

    block_model.FilmCondition(
        name='ConvectionBack',
        createStepName='heat_step',
        surface=region_back,
        filmCoeff=film_coeff,
        sinkTemperature=ambient_temp,
        definition=EMBEDDED_COEFF,
    )


# ------------------------------- GENERATING MESH --------------------------------------- #

    # SEED THE PART
    # minsizefactor is what is min element size compared to size
    # deviation size is how far is min deviation from curved surface
    mdb.models['Model-%d' %(q)].parts['Part-1'].seedPart(deviationFactor=0.1, 
        minSizeFactor=0.1, size=SEED_SIZE)
        
    # SET ELEMENT TYPE (matrix)
    # the three different elements types are fallbacks - if first doesnt work, try the next
    mdb.models['Model-%d' %(q)].parts['Part-1'].setElementType(elemTypes=(ElemType(
        elemCode=DC3D8, elemLibrary=STANDARD), ElemType(elemCode=DC3D6, 
        elemLibrary=STANDARD), ElemType(elemCode=DC3D4, elemLibrary=STANDARD)), 
        regions=(mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((MATRIX_POS,MATRIX_POS,MATRIX_POS), ), ), ))
        
    # SET ELEMENT TYPE (inclusions)
    # finds each inclusion based on their coords from centre coords array
    # then sets element type of each one
    for i in range (num_incl):
        mdb.models['Model-%d' %(q)].parts['Part-1'].setElementType(elemTypes=(ElemType(elemCode=DC3D8, elemLibrary=STANDARD), ElemType(elemCode=DC3D6, 
            elemLibrary=STANDARD), ElemType(elemCode=DC3D4, elemLibrary=STANDARD)), 
            regions=(mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i],y_coordinate[i]-0.2*RAD,z_coordinate[i]+0.2*RAD), ), ((x_coordinate[i],y_coordinate[i]+0.2*RAD,z_coordinate[i]+0.2*RAD), ), ), ))
            
    # SET MESH CONTROLS
    # tell abaqus to use tetrahedral elements with free meshing (on matrix)
    mdb.models['Model-%d' %(q)].parts['Part-1'].setMeshControls(elemShape=TET, regions=
        mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((MATRIX_POS, MATRIX_POS, MATRIX_POS), ), ), technique=FREE)
    
    # SET MESH CONTROLS
    # tell abaqus to mesh with tet elements using free meshing (in each inclusion)
    for i in range (num_incl):
        mdb.models['Model-%d' %(q)].parts['Part-1'].setMeshControls(elemShape=TET, regions=
    	    mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i],y_coordinate[i]-0.2*RAD,z_coordinate[i]+0.2*RAD), ), ((x_coordinate[i],y_coordinate[i]+0.2*RAD,z_coordinate[i]+0.2*RAD), ), ), technique=FREE)
        mdb.models['Model-%d' %(q)].parts['Part-1'].setMeshControls(elemShape=TET, regions=
    	    mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i],y_coordinate[i]-1.05*RAD,z_coordinate[i]), ), ((x_coordinate[i],y_coordinate[i]+1.05*RAD,z_coordinate[i]), ), ), technique=FREE)                

    # GENERATE MESH
    # literally just generate the mesh
    mdb.models['Model-%d' %(q)].parts['Part-1'].generateMesh()


 # ------ CREATE JOBS ----- #

    #LET'S CREATE JOBS 
    # mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    #     explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    #     memory=90, memoryUnits=PERCENTAGE, model='Model-%d' %(q), modelPrint=OFF, 
    #     multiprocessingMode=DEFAULT, name='Job-%d' %(q) , nodalOutputPrecision=SINGLE, 
    #     numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
    #     waitHours=0, waitMinutes=0)
    # mdb.jobs['Job-%d-%d' %(w,q)].writeInput()
    # mdb.jobs['Job-%d-%d' %(w,q) ].submit(consistencyChecking=OFF)    
    # mdb.jobs['Job-%d-%d' %(w,q) ].waitForCompletion()

    jobName   = 'HeatJob-%d' % q
    modelName = 'Model-%d' % q

    mdb.Job(name=jobName, model=modelName)
    mdb.jobs[jobName].submit()

# ---------- GEOMETRY CHECK ------------- #

for q in range (1,MAX_ITERATIONS+1):

    mdb.models['Model-%d' %(q)].parts['Part-1'].checkGeometry()

# prints out whether the partitioning worked as expected for each model



