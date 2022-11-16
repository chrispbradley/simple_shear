!> Main program
PROGRAM SimpleShearExample

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: PRESSURE_INTERPOLATION_TYPE=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  
  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: PRESSURE_BASIS_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: FIBRE_FIELD_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: MATERIAL_FIELD_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program types

  !Program variables
  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(CMISSIntg) :: numberOfGlobalXNodes,numberOfGlobalYNodes,numberOfGlobalZNodes
  INTEGER(CMISSIntg) :: numberOfDependentComponents,numberOfDimensions,numberOfGaussXi,numberOfMaterialsComponents, &
    & numberOfNodesPerElement,pressureComponent
  INTEGER(CMISSIntg) :: argumentLength,interpolationType,materialType,materialLaw,numberOfArguments,status
  INTEGER(CMISSIntg) :: dimensionIdx,nodeDomain,nodeNumber,xNodeIdx,zNodeIdx
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  LOGICAL :: directoryExists = .FALSE.
  LOGICAL :: usePressureBasis
  CHARACTER(LEN=255) :: commandArgument
 
  !OpenCMISS variables
  TYPE(cmfe_BasisType) :: basis,pressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: geometricField,fibreField,materialsField,dependentField,equationsSetField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver,linearSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !Generic OpenCMISS variables
  INTEGER(CMISSIntg) :: err

  !Setup problem 
  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 3) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalXElements
    IF(numberOfGlobalXElements<=0) CALL HandleError("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalYElements
    IF(numberOfGlobalYElements<=0) CALL HandleError("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalZElements
    IF(numberOfGlobalZElements<0) CALL HandleError("Invalid number of Z elements.")
    IF(numberOfArguments>=4) THEN
      CALL GET_COMMAND_ARGUMENT(4,commandArgument,argumentLength,status)
      IF(status>0) CALL HandleError("Error for command argument 4.")
      READ(commandArgument(1:argumentLength),*) interpolationType
      IF(interpolationType<=0) CALL HandleError("Invalid interpolation specification.")
      IF(numberOfArguments>=5) THEN
        CALL GET_COMMAND_ARGUMENT(5,commandArgument,argumentLength,status)
        IF(status>0) CALL HandleError("Error for command argument 5.")
        READ(commandArgument(1:argumentLength),*) materialType
        IF(materialType<=0) CALL HandleError("Invalid material type specification.")
      ELSE
        materialType=1
      ENDIF
    ELSE
      interpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
      materialType=1
    ENDIF
  ELSE
    !If there are not enough arguments default the problem specification
    numberOfGlobalXElements=2
    numberOfGlobalYElements=2
    numberOfGlobalZElements=2
    interpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    !interpolationType=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
    !interpolationType=CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION
    materialType=1
  ENDIF

  SELECT CASE(interpolationType)
  CASE(1,4,7)
    numberOfNodesPerElement=2
  CASE(2,8)
    numberOfNodesPerElement=3
  CASE(3,9)
    numberOfNodesPerElement=4
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")    
  END SELECT

  SELECT CASE(materialType)
  CASE(1)
    materialLaw=CMFE_EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE
    numberOfMaterialsComponents=3
    usePressureBasis=.FALSE.
  CASE(2)
    materialLaw=CMFE_EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE
    numberOfMaterialsComponents=2
    usePressureBasis=.TRUE.
  CASE(3)
    materialLaw=CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE
    numberOfMaterialsComponents=2
    usePressureBasis=.TRUE.
  CASE DEFAULT
    CALL HandleError("Invalid material type.")    
  END SELECT
  
  IF(numberOfGlobalZElements==0) THEN
    numberOfDimensions=2
    numberOfGaussXi=2
  ELSE
    numberOfDimensions=3
    numberOfGaussXi=3
  ENDIF
  numberOfGlobalXNodes=1+numberOfGlobalXElements*(numberOfNodesPerElement-1)
  numberOfGlobalYNodes=1+numberOfGlobalYElements*(numberOfNodesPerElement-1)
  numberOfGlobalZNodes=1+numberOfGlobalZElements*(numberOfNodesPerElement-1)
  IF(usePressureBasis) THEN
    numberOfDependentComponents=numberOfDimensions+1
    pressureComponent=numberOfDependentComponents
  ELSE
    numberOfDependentComponents=numberOfDimensions
    pressureComponent=0
  ENDIF

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Set all diganostic levels on for testing
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["FiniteElasticity_FiniteElementResidualEvaluate"],err)
  !Set output on
  CALL cmfe_OutputSetOn("SimpleShear",err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
  
  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)  

  !Create a rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions,err)
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"Region",err)
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL cmfe_Region_CreateFinish(region,err)

  !Define geometric basis
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  SELECT CASE(interpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_SIMPLEX_TYPE,err)
  END SELECT
  IF(numberOfDimensions==2) THEN
    CALL cmfe_Basis_NumberOfXiSet(basis,2,err)
    CALL cmfe_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ELSE
    CALL cmfe_Basis_NumberOfXiSet(basis,3,err)
    CALL cmfe_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ENDIF
  CALL cmfe_Basis_CreateFinish(basis,err)

  !Define pressure basis
  IF(usePressureBasis) THEN
    CALL cmfe_Basis_Initialise(pressureBasis,err)
    CALL cmfe_Basis_CreateStart(PRESSURE_BASIS_USER_NUMBER,context,pressureBasis,err)
    SELECT CASE(PRESSURE_INTERPOLATION_TYPE)
    CASE(1,2,3,4)
      CALL cmfe_Basis_TypeSet(pressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
    CASE(7,8,9)
      CALL cmfe_Basis_TypeSet(pressureBasis,CMFE_BASIS_SIMPLEX_TYPE,err)
    END SELECT
    IF(numberOfDimensions==2) THEN
      CALL cmfe_Basis_NumberOfXiSet(pressureBasis,2,err)
      CALL cmfe_Basis_InterpolationXiSet(pressureBasis,[PRESSURE_INTERPOLATION_TYPE,PRESSURE_INTERPOLATION_TYPE],err)
      IF(numberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(pressureBasis,[numberOfGaussXi,numberOfGaussXi],err)
      ENDIF
    ELSE
      CALL cmfe_Basis_NumberOfXiSet(pressureBasis,3,err)
      CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
        & [PRESSURE_INTERPOLATION_TYPE,PRESSURE_INTERPOLATION_TYPE,PRESSURE_INTERPOLATION_TYPE],err)
      IF(numberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(pressureBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
      ENDIF
    ENDIF
    CALL cmfe_Basis_CreateFinish(pressureBasis,err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  IF(usePressureBasis) THEN
    CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,[basis,pressureBasis],err)
  ELSE
    CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,[basis],err)
  ENDIF
  !Define the mesh on the region
  IF(numberOfDimensions==2) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  CALL cmfe_Field_VariableLabelSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL cmfe_Field_ScalingTypeSet(geometricField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,err)
  CALL cmfe_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(fibreField,err)
  CALL cmfe_Field_CreateStart(FIBRE_FIELD_USER_NUMBER,region,fibreField,err)
  CALL cmfe_Field_TypeSet(fibreField,CMFE_FIELD_FIBRE_TYPE,err)
  CALL cmfe_Field_DecompositionSet(fibreField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(fibreField,geometricField,err)
  CALL cmfe_Field_VariableLabelSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",err)
  CALL cmfe_Field_CreateFinish(fibreField,err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_Field_CreateStart(DEPENDENT_FIELD_USER_NUMBER,region,dependentField,err)
  CALL cmfe_Field_TypeSet(dependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,err)
  CALL cmfe_Field_DecompositionSet(dependentField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(dependentField,geometricField,err)
  CALL cmfe_Field_DependentTypeSet(dependentField,CMFE_FIELD_DEPENDENT_TYPE,err)
  CALL cmfe_Field_NumberOfVariablesSet(dependentField,2,err)
  CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfDependentComponents,err)
  IF(usePressureBasis) THEN
    !Set the pressure to be nodally based and use the second mesh component if required
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,pressureComponent, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,pressureComponent, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,pressureComponent,2,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,pressureComponent,2,err)
  END IF
  CALL cmfe_Field_CreateFinish(dependentField,err)

  !Create the material field
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_Field_CreateStart(MATERIAL_FIELD_USER_NUMBER,region,materialsField,err)
  CALL cmfe_Field_TypeSet(materialsField,CMFE_FIELD_MATERIAL_TYPE,err)
  CALL cmfe_Field_DecompositionSet(materialsField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(materialsField,geometricField,err)
  CALL cmfe_Field_NumberOfVariablesSet(materialsField,1,err)
  CALL cmfe_Field_VariableLabelSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL cmfe_Field_NumberOfComponentsSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents,err)
  CALL cmfe_Field_CreateFinish(materialsField,err)

  !Create the equations_set
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,materialLaw],EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !Create the equations set dependent field
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the equations set material field 
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIAL_FIELD_USER_NUMBER,materialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Set Mooney-Rivlin constants c10 and c01 to 0.5 and 0.0 respectively. Third value is kappa (bulk modulus ???)
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,0.0_CMISSRP,err)
  IF(numberOfMaterialsComponents==3) THEN
    CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 3,10000.0_CMISSRP,err)
  ENDIF

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !Initialise dependent field from undeformed geometry and displacement bcs
  DO dimensionIdx=1,numberOfDimensions
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & dimensionIdx,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,dimensionIdx,err)
  ENDDO !dimensionIdx

  !Define the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
  CALL cmfe_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  CALL cmfe_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  !Set x=LENGTH nodes to alpha% x-displacement, no displacement in y- and z-direction
  DO zNodeIdx=1,numberOfGlobalZNodes
    DO xNodeIdx=1,numberOfGlobalXNodes
      !Fix the bottom nodes in all directions
      nodeNumber=xNodeIdx+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !x-direction
        CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
        !y-direction
        CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
          & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
        IF(numberOfDimensions==3) THEN
          !z-direction
          CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
            & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
        ENDIF
      ENDIF
      !Fix the top nodes to 10% x-displacement and fixing the other directions
      nodeNumber=xNodeIdx+numberOfGlobalXNodes*(numberOfGlobalYNodes-1)+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !x-direction
        CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & CMFE_BOUNDARY_CONDITION_FIXED,0.1_CMISSRP*WIDTH,err)
        !y-direction
        CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
          & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
        IF(numberOfDimensions==3) THEN
          !z-direction
          CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
            & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
        ENDIF
      ENDIF
    ENDDO !xNodeIdx
  ENDDO !zNodeIdx
  
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL cmfe_Problem_Solve(problem,err)

  INQUIRE(FILE="./results",EXIST=directoryExists)
  IF (.NOT.directoryExists) THEN
    CALL EXECUTE_COMMAND_LINE("mkdir ./results")
  ENDIF

  !Output solution
  CALL cmfe_Fields_Initialise(fields,err)
  CALL cmfe_Fields_Create(region,fields,err)
  CALL cmfe_Fields_NodesExport(fields,"./results/SimpleShear","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(fields,"./results/SimpleShear","FORTRAN",err)
  CALL cmfe_Fields_Finalise(fields,err)

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError

END PROGRAM SimpleShearExample

