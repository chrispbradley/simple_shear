!> Main program
PROGRAM SimpleShearExample

  USE OpenCMISS

  IMPLICIT NONE
  
  !Test program parameters

  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=1.0_OC_RP
  
  INTEGER(OC_Intg), PARAMETER :: PRESSURE_INTERPOLATION_TYPE=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  
  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: PRESSURE_BASIS_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: FIBRE_FIELD_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: MATERIAL_FIELD_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program types

  !Program variables
  INTEGER(OC_Intg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(OC_Intg) :: numberOfGlobalXNodes,numberOfGlobalYNodes,numberOfGlobalZNodes
  INTEGER(OC_Intg) :: numberOfDependentComponents,numberOfDimensions,numberOfGaussXi,numberOfMaterialsComponents, &
    & numberOfNodesPerElement,pressureComponent
  INTEGER(OC_Intg) :: argumentLength,interpolationType,materialType,materialLaw,numberOfArguments,status
  INTEGER(OC_Intg) :: dimensionIdx,nodeDomain,nodeNumber,xNodeIdx,zNodeIdx
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex
  INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber
  LOGICAL :: directoryExists = .FALSE.
  LOGICAL :: usePressureBasis
  CHARACTER(LEN=255) :: commandArgument
 
  !OpenCMISS variables
  TYPE(OC_BasisType) :: basis,pressureBasis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: geometricField,fibreField,materialsField,dependentField,equationsSetField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: solver,linearSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !Generic OpenCMISS variables
  INTEGER(OC_Intg) :: err

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
      interpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
      materialType=1
    ENDIF
  ELSE
    !If there are not enough arguments default the problem specification
    numberOfGlobalXElements=2
    numberOfGlobalYElements=2
    numberOfGlobalZElements=2
    interpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    !interpolationType=OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
    !interpolationType=OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION
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
    materialLaw=OC_EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE
    numberOfMaterialsComponents=3
    usePressureBasis=.FALSE.
  CASE(2)
    materialLaw=OC_EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE
    numberOfMaterialsComponents=2
    usePressureBasis=.TRUE.
  CASE(3)
    materialLaw=OC_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE
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
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Set all diganostic levels on for testing
  CALL OC_DiagnosticsSetOn(OC_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["FiniteElasticity_FiniteElementResidualEvaluate"],err)
  !Set output on
  CALL OC_OutputSetOn("SimpleShear",err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)
  
  !Get the number of computational nodes and this computational node number
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)  

  !Create a rectangular cartesian coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions,err)
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the coordinate system to the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"Region",err)
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL OC_Region_CreateFinish(region,err)

  !Define geometric basis
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  SELECT CASE(interpolationType)
  CASE(1,2,3,4)
    CALL OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL OC_Basis_TypeSet(basis,OC_BASIS_SIMPLEX_TYPE,err)
  END SELECT
  IF(numberOfDimensions==2) THEN
    CALL OC_Basis_NumberOfXiSet(basis,2,err)
    CALL OC_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ELSE
    CALL OC_Basis_NumberOfXiSet(basis,3,err)
    CALL OC_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ENDIF
  CALL OC_Basis_CreateFinish(basis,err)

  !Define pressure basis
  IF(usePressureBasis) THEN
    CALL OC_Basis_Initialise(pressureBasis,err)
    CALL OC_Basis_CreateStart(PRESSURE_BASIS_USER_NUMBER,context,pressureBasis,err)
    SELECT CASE(PRESSURE_INTERPOLATION_TYPE)
    CASE(1,2,3,4)
      CALL OC_Basis_TypeSet(pressureBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
    CASE(7,8,9)
      CALL OC_Basis_TypeSet(pressureBasis,OC_BASIS_SIMPLEX_TYPE,err)
    END SELECT
    IF(numberOfDimensions==2) THEN
      CALL OC_Basis_NumberOfXiSet(pressureBasis,2,err)
      CALL OC_Basis_InterpolationXiSet(pressureBasis,[PRESSURE_INTERPOLATION_TYPE,PRESSURE_INTERPOLATION_TYPE],err)
      IF(numberOfGaussXi>0) THEN
        CALL OC_Basis_QuadratureNumberOfGaussXiSet(pressureBasis,[numberOfGaussXi,numberOfGaussXi],err)
      ENDIF
    ELSE
      CALL OC_Basis_NumberOfXiSet(pressureBasis,3,err)
      CALL OC_Basis_InterpolationXiSet(pressureBasis, &
        & [PRESSURE_INTERPOLATION_TYPE,PRESSURE_INTERPOLATION_TYPE,PRESSURE_INTERPOLATION_TYPE],err)
      IF(numberOfGaussXi>0) THEN
        CALL OC_Basis_QuadratureNumberOfGaussXiSet(pressureBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
      ENDIF
    ENDIF
    CALL OC_Basis_CreateFinish(pressureBasis,err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  IF(usePressureBasis) THEN
    CALL OC_GeneratedMesh_BasisSet(generatedMesh,[basis,pressureBasis],err)
  ELSE
    CALL OC_GeneratedMesh_BasisSet(generatedMesh,[basis],err)
  ENDIF
  !Define the mesh on the region
  IF(numberOfDimensions==2) THEN
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL OC_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (default is geometry)
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  CALL OC_Field_VariableLabelSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL OC_Field_ScalingTypeSet(geometricField,OC_FIELD_ARITHMETIC_MEAN_SCALING,err)
  CALL OC_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create a fibre field and attach it to the geometric field
  CALL OC_Field_Initialise(fibreField,err)
  CALL OC_Field_CreateStart(FIBRE_FIELD_USER_NUMBER,region,fibreField,err)
  CALL OC_Field_TypeSet(fibreField,OC_FIELD_FIBRE_TYPE,err)
  CALL OC_Field_DecompositionSet(fibreField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(fibreField,geometricField,err)
  CALL OC_Field_VariableLabelSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,"Fibre",err)
  CALL OC_Field_CreateFinish(fibreField,err)

  !Create the dependent field
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_Field_CreateStart(DEPENDENT_FIELD_USER_NUMBER,region,dependentField,err)
  CALL OC_Field_TypeSet(dependentField,OC_FIELD_GEOMETRIC_GENERAL_TYPE,err)
  CALL OC_Field_DecompositionSet(dependentField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(dependentField,geometricField,err)
  CALL OC_Field_DependentTypeSet(dependentField,OC_FIELD_DEPENDENT_TYPE,err)
  CALL OC_Field_NumberOfVariablesSet(dependentField,2,err)
  CALL OC_Field_VariableTypesSet(dependentField,[OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_T_VARIABLE_TYPE],err)
  CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  CALL OC_Field_NumberOfComponentsSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,numberOfDependentComponents,err)
  IF(usePressureBasis) THEN
    !Set the pressure to be nodally based and use the second mesh component if required
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,pressureComponent, &
      & OC_FIELD_NODE_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,pressureComponent, &
      & OC_FIELD_NODE_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,pressureComponent,2,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,pressureComponent,2,err)
  END IF
  CALL OC_Field_CreateFinish(dependentField,err)

  !Create the material field
  CALL OC_Field_Initialise(materialsField,err)
  CALL OC_Field_CreateStart(MATERIAL_FIELD_USER_NUMBER,region,materialsField,err)
  CALL OC_Field_TypeSet(materialsField,OC_FIELD_MATERIAL_TYPE,err)
  CALL OC_Field_DecompositionSet(materialsField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(materialsField,geometricField,err)
  CALL OC_Field_NumberOfVariablesSet(materialsField,1,err)
  CALL OC_Field_VariableLabelSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL OC_Field_NumberOfComponentsSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents,err)
  CALL OC_Field_CreateFinish(materialsField,err)

  !Create the equations_set
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,materialLaw],EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  !Create the equations set dependent field
  CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the equations set material field 
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIAL_FIELD_USER_NUMBER,materialsField,err)
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Set Mooney-Rivlin constants c10 and c01 to 0.5 and 0.0 respectively. Third value is kappa (bulk modulus ???)
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,0.5_OC_RP,err)
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,0.0_OC_RP,err)
  IF(numberOfMaterialsComponents==3) THEN
    CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & 3,10000.0_OC_RP,err)
  ENDIF

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !Initialise dependent field from undeformed geometry and displacement bcs
  DO dimensionIdx=1,numberOfDimensions
    CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & dimensionIdx,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,dimensionIdx,err)
  ENDDO !dimensionIdx

  !Define the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & OC_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
  CALL OC_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolversCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT,err)
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(solver,OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  CALL OC_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  !Set x=LENGTH nodes to alpha% x-displacement, no displacement in y- and z-direction
  DO zNodeIdx=1,numberOfGlobalZNodes
    DO xNodeIdx=1,numberOfGlobalXNodes
      !Fix the bottom nodes in all directions
      nodeNumber=xNodeIdx+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !x-direction
        CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
        !y-direction
        CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
          & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
        IF(numberOfDimensions==3) THEN
          !z-direction
          CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
            & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
        ENDIF
      ENDIF
      !Fix the top nodes to 10% x-displacement and fixing the other directions
      nodeNumber=xNodeIdx+numberOfGlobalXNodes*(numberOfGlobalYNodes-1)+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
      CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
      IF(nodeDomain==computationalNodeNumber) THEN
        !x-direction
        CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
          & OC_BOUNDARY_CONDITION_FIXED,0.1_OC_RP*WIDTH,err)
        !y-direction
        CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,2, &
          & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
        IF(numberOfDimensions==3) THEN
          !z-direction
          CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,3, &
            & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
        ENDIF
      ENDIF
    ENDDO !xNodeIdx
  ENDDO !zNodeIdx
  
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL OC_Problem_Solve(problem,err)

  INQUIRE(FILE="./results",EXIST=directoryExists)
  IF (.NOT.directoryExists) THEN
    CALL EXECUTE_COMMAND_LINE("mkdir ./results")
  ENDIF

  !Output solution
  CALL OC_Fields_Initialise(fields,err)
  CALL OC_Fields_Create(region,fields,err)
  CALL OC_Fields_NodesExport(fields,"./results/SimpleShear","FORTRAN",err)
  CALL OC_Fields_ElementsExport(fields,"./results/SimpleShear","FORTRAN",err)
  CALL OC_Fields_Finalise(fields,err)

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError

END PROGRAM SimpleShearExample

