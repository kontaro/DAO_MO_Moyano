package dao_singobj_ls;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import IMRT_Base.*;
import gurobi.GRBException;

/**
*
* @author Denisse and Mauricio
*/

public class DAO_SingObj_LS {
	
	public static int numOrgans;
    public static int numAngles;                   //#Beams in TreatmentPlan
    public static Organs[] o;
    public static String pathFile = "";
    public static int option;
    public static int selectionCriterion;
    public static int[][] Vx; 				       //Vx index. % of the volume receiving more than x% of the prescribed dose
    						  				       //One row per organ (usually only for target regions)
    public static int[][] Dx; 				       //Dx index. Minimum dose of the hottest x% of the volume.
    						  				       //One row per organ (usually only for OAR regions)
    public static boolean randomInitialAngles;
    public static int[] initialBACs;
    public static int[][] beamletSet;
    
    public static int num_apertures = 5;            // # DAO : Number of apertures per station
	public static int init_intensity;               // # DAO : Initial intensity for apertures
    public static int max_intensity;                // # DAO : Maximum intensity for apertures
    public static int max_delta;                    // # DAO : Delta máximo para la variación de intesidad por apertura
    public static int max_iter;              		// # DAO : Cantidad máxima de iteraciones del algoritmo
    public static int max_time;                     // # DAO : tiempo máximo de ejecución del algoritmo
    public static int seed;                         // # DAO : semilla para partir la busqueda local
    public static int step_intensity;				// # DAO : Step size for aperture intensity (2)
    public static Vector<int[][]> initial_aperture_shape;  // # DAO : Conjunto de formas de las aperturas según restricciones (rango sup y rango inf)
    public static Vector<int[][]> apertureShapes;   // # DAO: Aperturas segun restricciones del problema (en forma de matriz)
    public static int [] zmax,zmin,dd;
    public static int totalbmlt;
    
	public static void main(String[] args) throws IOException{
		
		// Lectura de archivos y creacion de DDM  y seteo de parametors
		
		String inputFile = args[0];
		
	    readInputFile("DAO_SingObj_LS/"+inputFile);
	    beamletSet = new int[numAngles][4];
	    
	    TreatmentPlan solution = new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
	    int[] selAngles = new int[]{70,65,67,71,63};
	    int totalbmlt =336;
	    //Constructor DAO for DDM
	    //long time=System.currentTimeMillis();
	    String beamInfoDir = pathFile+"beamsInfo.txt";
	    
	    DDM M; M = new DDM(o, initialBACs, pathFile,selAngles,totalbmlt);
	    
	    
	    
	    beamletSet = M.loadNumBixels(beamInfoDir);
	    //System.out.println(System.currentTimeMillis()-time);
	    zmin = new int[]{70,65,65};
	    zmax = new int[]{1000,55,50};
	    double[] weight = new double[]{0.5,0.35,0.15};
	    dd = new int[]{70,55,50};
	    double[] w = new double[]{0.5,0.35,0.15}; 
	    solution.setWeights(w);
	    setBeams(solution, pathFile);
	    setInitialApertureShape(solution);
	    generateFirstSolution(solution);
	   
	    
	    
	    //Algoritmo a ejecutar
	    
	    
	    switch(selectionCriterion){
		    case 0://Local Search con modelo Matematico
				try {
					
					LSmath(M,1,selectionCriterion);
				} catch (GRBException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		    	break;
		    case 1:
				try {
					solution=Two_Phase(M,2,selectionCriterion);
				} catch (GRBException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
		    	break;
		    case 2:
		    	MO(M);
		    	break;
		    default:
		    	
	    	
	    }

	}
	
	
	//-------------------------------------Algoritmos que utilizan 
	
	/**
	 * Algoritmo lineal el cual genera una solucion inicial y realiza una busqueda lcoal, para posteriormente usar un solver para encontrar las 
	 * mejores intesidades
	 * @param M
	 * @param initial
	 * @param moviment
	 * @return
	 * @throws IOException
	 * @throws GRBException
	 */
	public static TreatmentPlan Two_Phase(DDM M,int initial,int moviment) throws IOException, GRBException {
	    TreatmentPlan solution = new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
	    //Generacion random
	    double[] w = new double[]{0.5,0.35,0.15}; 
	    int[] selAngles = new int[]{70,65,67,71,63};
	    solution.setWeights(w);
	    setBeams(solution, pathFile);
	    setInitialApertureShape(solution);
	    
	    switch(initial) {
	    	case 1:
	    		generateFirstSolution(solution);
	    		break;
	    	case 2:
	    		generateFirstSolutionWithFixedIntensity(solution);
	    		break;
	    	
	    	case 3:
	    		generateFirstSolutionOpen(solution);
	    		break;
	    	default:
	    		generateFirstSolutionOpenWithFixedIntensity(solution);
	    		break;
	    		
	    }
	    
	    solution.aperturesToIntensityMatrixPerStationwithoutIntensty();
	    solution.ApertureMatrixStationToIntensityMatrixPlan();
	    solution.IntensityMatrixPlan();

	    evaluateSolution(solution,M,o,zmin,zmax,solution.weights);

	    // Local search
	    solution.updateSolDAO(LocalSearch(solution,M,1,5,moviment));
	    
	    //Solver
	    Gurobi_Solver newModel=new Gurobi_Solver(solution,M,o,selAngles,dd,w);
	    
	    solution.setIntensity(newModel.newIntensity);
	    solution.IntensityMatrixPlan();
	    evaluateSolution(solution,M,o,zmin,zmax,solution.weights);
	    System.out.println(solution.singleObjectiveValue);
	    
	    
	    
	    return solution;
	}
	
	public static void MO(DDM M) throws IOException {
		ArrayList<TreatmentPlan> notDominated=new ArrayList<TreatmentPlan>();
		TreatmentPlan solution = new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		double[] w = new double[]{0.5,0.35,0.15}; 
	    int[] selAngles = new int[]{70,65,67,71,63};
	    solution.setWeights(w);
	    setBeams(solution, pathFile);
	    setInitialApertureShape(solution);
		generateFirstSolutionWithFixedIntensity(solution);
		
		notDominated=ParetoLocalSearch(solution,M);
		System.out.println(notDominated.size());
		printFront(notDominated);
	}
	
	

	
	public static TreatmentPlan LSmath(DDM M,int initial,int moviment) throws IOException, GRBException {
	    TreatmentPlan solution = new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);

	    double[] w = new double[]{0.5,0.35,0.15}; 
	    int[] selAngles = new int[]{70,65,67,71,63};
	    solution.setWeights(w);
	    setBeams(solution, pathFile);
	    setInitialApertureShape(solution);
	    moviment=1;
	    switch(initial) {
	    	case 1:
	    		generateFirstSolutionWithFixedIntensity(solution);
	    		break;
	    	default:
	    		generateFirstSolutionOpenWithFixedIntensity(solution);
	    		break;
	    		
	    }
	   
	    solution.aperturesToIntensityMatrixPerStationwithoutIntensty();
	    
	    
	    solution.ApertureMatrixStationToIntensityMatrixPlan();
	    solution.IntensityMatrixPlan();
	    Gurobi_Solver newModel=new Gurobi_Solver(solution,M,o,selAngles,dd,w);
	    solution.setIntensity(newModel.newIntensity);
	    solution.setOF(newModel.objVal);
	    System.out.print("initial value:"+solution.singleObjectiveValue+" ");
//	    printapertureBeamlet(solution);
	    // Local search
	    solution.updateSolDAO(LocalSearchWithSolver(solution,M,newModel,3,5, moviment));

	    System.out.println("final value:"+solution.singleObjectiveValue);
	    
	    
	    
	    return solution;
	}

	
	//-----------------------Algoritmos monoObjetivo------------------------------------------
	
	/**
	 * Va generando soluciones random y evalua si es mejor que la actual
	 * @param sol
	 * @param M
	 * @param numberGenerations
	 * @return
	 */
	public static TreatmentPlan RandomSearch(TreatmentPlan sol,DDM M,int numberGenerations) {
    	TreatmentPlan bestTP=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
    	bestTP.updateSolDAO(sol);

	    
	    for(int i=0;i<numberGenerations;i++) {
	    	TreatmentPlan actualTP=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			actualTP.updateSolDAO(sol);
			generateFirstSolution(actualTP);
			actualTP.aperturesToIntensityMatrixPerStation();
		    actualTP.intensityMatrixStationToIntensityMatrixPlan();
		    evaluateSolution(actualTP,M,o,zmin,zmax,actualTP.weights);
		    if(actualTP.singleObjectiveValue<bestTP.singleObjectiveValue) {
		    	bestTP.updateSolDAO(actualTP);
		    }
	    }
	   
		
		return bestTP;
	}
	

	

	
	//-----------------------Algoritmos basados en busqueda local--------------------------
	
	public static TreatmentPlan LocalSearch(TreatmentPlan sol,DDM M,int conditionCase,float conditionValue,int moviment) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood;
		int mov=1;
		int iter=0;
		float iterMax=conditionValue;
		boolean condition=true;
		boolean update=false;
		float maxTime=conditionValue;
		float time=System.currentTimeMillis();
		/**ITERACION
		 * Se itera hasta:
		 * Convergencia
		 * Tiempo maximo 
		 * Iteraciones maxima - Actualmente se usa esta
		 */
		
		switch(conditionCase) {
			case 1:
				while(condition) {
					//System.out.println("entre");
					if(moviment ==1) {
						neigborhood=bestMove(BestSolution);
					}else {
						neigborhood=firstMove(BestSolution,M);
					}
					if(neigborhood.size()>0) {
						actualSolution.updateSolDAO(bestNeigbor(neigborhood,mov,M));
						//System.out.println(actualSolution.singleObjectiveValue);
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
						}
						iter=iter+1;
						
						if(iterMax<iter) {
							condition=false;
						}
					}
					else {
						condition=false;
					}
				}
				break;
			case 2:
				while(condition) {
					if(moviment ==1) {
						neigborhood=bestMove(BestSolution);
					}else {
						neigborhood=firstMove(BestSolution,M);
					}
					if(neigborhood.size()>0) {
						actualSolution.updateSolDAO(bestNeigbor(neigborhood,mov,M));
						//System.out.println(actualSolution.singleObjectiveValue);
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
						}

						if((System.currentTimeMillis()-time)>maxTime) {
							condition=false;
						}
					}
					else {
						condition=false;
					}
					
				}
				break;
			case 3:
				while(condition) {
					update=false;
					neigborhood=bestMove(BestSolution);
					actualSolution.updateSolDAO(bestNeigbor(neigborhood,mov,M));
					//System.out.println(actualSolution.singleObjectiveValue);
					if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
						BestSolution.updateSolDAO(actualSolution);
						update=true;
					}
					if(!update) {
						condition=false;
					}

				}
				break;
			default:
				while(condition) {
					
					neigborhood=bestMove(BestSolution);
					actualSolution.updateSolDAO(bestNeigbor(neigborhood,mov,M));
					//System.out.println(actualSolution.singleObjectiveValue);
					if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
						BestSolution.updateSolDAO(actualSolution);
					}
					iter=iter+1;
					
					if(iterMax<iter) {
						condition=false;
					}
				}
				break;
				
		
		}
		
	
		
		return BestSolution;
	}
	
	/**
	 * Es una busqueda local que al generar un vecino lo prueba lo evalua con el solver
	 * @param sol
	 * @param M
	 * @param conditionCase
	 * @param conditionValue
	 * @param moviment
	 * @return
	 * @throws IOException
	 */
	public static TreatmentPlan LocalSearchWithSolver(TreatmentPlan sol,DDM M,Gurobi_Solver GS,int conditionCase,float conditionValue,int moviment) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood;
		int mov=1;
		int iter=0;
		float iterMax=conditionValue;
		boolean condition=true;
		boolean update=false;
		float maxTime=conditionValue;
		float time=System.currentTimeMillis();
		/**ITERACION
		 * Se itera hasta:
		 * Convergencia
		 * Tiempo maximo 
		 * Iteraciones maxima - Actualmente se usa esta
		 */
		
		switch(conditionCase) {
			case 1:
				while(condition) {
					
					if(moviment ==1) {
						neigborhood=bestMove(BestSolution);
					}else {
						neigborhood=firstMoveSolver(BestSolution,M,GS);
					}
					if(neigborhood.size()>0) {
						actualSolution.updateSolDAO(bestNeigbor(neigborhood,mov,M));
						//System.out.println(actualSolution.singleObjectiveValue);
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
						}
						iter=iter+1;
						
						if(iterMax<iter) {
							condition=false;
						}
					}
					else {
						condition=false;
					}
				}
				break;
			case 2:
				while(condition) {
					if(moviment ==1) {
						neigborhood=bestMove(BestSolution);
					}else {
						neigborhood=firstMoveSolver(BestSolution,M,GS);
					}
					if(neigborhood.size()>0) {
						actualSolution.updateSolDAO(bestNeigbor(neigborhood,mov,M));
						//System.out.println(actualSolution.singleObjectiveValue);
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
						}

						if((System.currentTimeMillis()-time)>maxTime) {
							condition=false;
						}
					}
					else {
						condition=false;
					}
					
				}
				break;
			case 3:
				while(condition) {
					update=false;
					if(moviment ==1) {
						neigborhood=bestMove(BestSolution);
					}else {
						neigborhood=firstMoveSolver(BestSolution,M,GS);
					}
					if(neigborhood.size()>0) {
						if(moviment==1) {
							actualSolution.updateSolDAO(bestNeigborSolver(neigborhood,mov,M,GS));
							//System.out.println(actualSolution.singleObjectiveValue);
							if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
								BestSolution.updateSolDAO(actualSolution);
								update=true;
								System.out.println(BestSolution.singleObjectiveValue);
							}
							if(!update) {
								condition=false;
							}
						}else {
							actualSolution.updateSolDAO(bestNeigbor(neigborhood,mov,M));
							//System.out.println(actualSolution.singleObjectiveValue);
							if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
								BestSolution.updateSolDAO(actualSolution);
								update=true;
								System.out.println(BestSolution.singleObjectiveValue);
							}
							if(!update) {
								condition=false;
							}
						}
					}else {
						condition=false;
					}
					

				}
				break;
			default:
				while(condition) {
					
					neigborhood=bestMove(BestSolution);
					actualSolution.updateSolDAO(bestNeigbor(neigborhood,mov,M));
					//System.out.println(actualSolution.singleObjectiveValue);
					if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
						BestSolution.updateSolDAO(actualSolution);
					}
					iter=iter+1;
					
					if(iterMax<iter) {
						condition=false;
					}
				}
				break;
				
		
		}
		
	
		
		return BestSolution;
	}
	
	
	//-----------------------Algoritmos Multi Objetivos basado en busqueda local-----------
	public static ArrayList<TreatmentPlan> ParetoLocalSearch(TreatmentPlan initialSolution,DDM M) throws IOException {
		ArrayList<TreatmentPlan> Front=new ArrayList<TreatmentPlan>();
		ArrayList<TreatmentPlan> neigborhood=new ArrayList<TreatmentPlan>();
		ArrayList<TreatmentPlan> visited=new ArrayList<TreatmentPlan>();
		int[] selAngles = new int[]{70,65,67,71,63};
		boolean allVisited=false;
		// Generar solucion Inicial
		Front.add(initialSolution);
		System.out.println(Front.size());
		//generar vecindad de los nodos no visitados - First Phase
		while(allVisited!=true) {
			
			
			for(TreatmentPlan nodo:Front){
				if(!visited.contains(nodo)) {
					neigborhood.addAll(bestMove(nodo));
					visited.add(nodo);
				}
				
			}
			System.out.println(neigborhood.size());
			updateFrontier(neigborhood,Front,M);
			
			for(TreatmentPlan nodo:Front){
				if(!visited.contains(nodo)) {
					 try {
						Gurobi_Solver newModel=new Gurobi_Solver(nodo,M,o,selAngles,dd,nodo.weights);
						nodo.setIntensity(newModel.newIntensity);
						nodo.IntensityMatrixPlan();
					    evaluateSolution(nodo,M,o,zmin,zmax,nodo.weights);
					} catch (GRBException e) {
						// TODO Auto-generated catch block
						
					}
					
				}
				
			}
			allVisited=true;
			if(visited.containsAll(Front))
                allVisited=true;
		}
		
		return Front;
	}
	//_----------------------Movimientos para busqueda local-------------------------------
	/* Estas funciones toman una solucion, y utilizando un movimiento de vencidad , 
	 * modifican la solucion a la de un vecino.
	 *  
	 */
	
	/**
	 * Esta funcion busca agrandar o achicar la apertura de una hoja, genera todos los vecinos posibles para luego ser evaluados 
	 * La apertura se escoge de forma random
	 * La hoja se escoge de forma random
	 * @return
	 * @throws IOException 
	 */
	public static ArrayList<TreatmentPlan> bestMove(TreatmentPlan sol) throws IOException {
		ArrayList<TreatmentPlan> neighborhood=new ArrayList<TreatmentPlan>();
		Random r = new Random();
		int angleSize=sol.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize);
		int apertureChoose= r.nextInt(sol.apertures.get(angleChoose).size()-1);
		List<Aperture> listApertures = sol.apertures.elementAt(angleSize);

		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		for(int i=0;i<sizeAperture;i++){
			//int[][] aperture = listApertures.get(apertureChoose).aperture;
			TreatmentPlan neigbor=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor.updateSolDAO(sol);
			neigbor.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 0, true);

			
			TreatmentPlan neigbor1=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor1.updateSolDAO(sol);
			neigbor1.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 1, true);
			
			TreatmentPlan neigbor2=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor2.updateSolDAO(sol);
			neigbor2.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 0, false);
			
			TreatmentPlan neigbor3=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor3.updateSolDAO(sol);
			neigbor3.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 1, false);
			//neigbor.changeAperture(aperture, angleChoose, apertureChoose);
			
			if(isFeasible(neigbor)){neighborhood.add(neigbor);}
			if(isFeasible(neigbor1)){neighborhood.add(neigbor1);}
			if(isFeasible(neigbor2)){neighborhood.add(neigbor2);}
			if(isFeasible(neigbor3)){neighborhood.add(neigbor3);}	
		}

		return neighborhood;
	}
	public static ArrayList<TreatmentPlan> firstMove(TreatmentPlan sol,DDM M) throws IOException {
		ArrayList<TreatmentPlan> neighborhood=new ArrayList<TreatmentPlan>();
		Random r = new Random();
		int angleSize=sol.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize);
		int apertureChoose= r.nextInt(sol.apertures.get(angleChoose).size()-1);
		List<Aperture> listApertures = sol.apertures.elementAt(angleSize);

		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		for(int i=0;i<sizeAperture;i++){
			
			TreatmentPlan neigbor=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor.updateSolDAO(sol);
			neigbor.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 0, true);
			if(isFeasible(neigbor)){
				neigbor.aperturesToIntensityMatrixPerStation();
				neigbor.intensityMatrixStationToIntensityMatrixPlan();
				evaluateSolution(neigbor,M,o,zmin,zmax,neigbor.weights);
				if(neigbor.singleObjectiveValue<sol.singleObjectiveValue) {
					neighborhood.add(neigbor);
					return neighborhood;
				}
			}
			
			TreatmentPlan neigbor1=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor1.updateSolDAO(sol);
			neigbor1.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 1, true);
			if(isFeasible(neigbor1)){
				neigbor1.aperturesToIntensityMatrixPerStation();
				neigbor1.intensityMatrixStationToIntensityMatrixPlan();
				evaluateSolution(neigbor1,M,o,zmin,zmax,neigbor1.weights);
				if(neigbor1.singleObjectiveValue<sol.singleObjectiveValue) {
					neighborhood.add(neigbor1);
					return neighborhood;
				}
			}
			TreatmentPlan neigbor2=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor2.updateSolDAO(sol);
			neigbor2.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 0, false);
			if(isFeasible(neigbor2)){
				neigbor2.aperturesToIntensityMatrixPerStation();
				neigbor2.intensityMatrixStationToIntensityMatrixPlan();
				evaluateSolution(neigbor2,M,o,zmin,zmax,neigbor2.weights);
				if(neigbor2.singleObjectiveValue<sol.singleObjectiveValue) {
					neighborhood.add(neigbor2);
					return neighborhood;
				}
			}
			TreatmentPlan neigbor3=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor3.updateSolDAO(sol);
			neigbor3.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 1, false);
			//neigbor.changeAperture(aperture, angleChoose, apertureChoose);
			if(isFeasible(neigbor3)){
				neigbor3.aperturesToIntensityMatrixPerStation();
				neigbor3.intensityMatrixStationToIntensityMatrixPlan();
				evaluateSolution(neigbor3,M,o,zmin,zmax,neigbor3.weights);
				if(neigbor3.singleObjectiveValue<sol.singleObjectiveValue) {
					neighborhood.add(neigbor3);
					return neighborhood;
				}
			}
	
		}

		return neighborhood;
	}
	
	/**
	 * First Move 
	 * @param sol
	 * @param M
	 * @return El primer vecino que sea mejor que la solucion actual
	 * @throws IOException
	 */
	public static ArrayList<TreatmentPlan> firstMoveSolver(TreatmentPlan sol,DDM M,Gurobi_Solver GS) throws IOException {
		ArrayList<TreatmentPlan> neighborhood=new ArrayList<TreatmentPlan>();
		Random r = new Random();
		int angleSize=sol.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize);
		int apertureChoose= r.nextInt(sol.apertures.get(angleChoose).size()-1);
		List<Aperture> listApertures = sol.apertures.elementAt(angleSize);
		//System.out.println("entre");
		int count=0;
		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		for(int i=0;i<sizeAperture;i++){
			
			TreatmentPlan neigbor=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor.updateSolDAO(sol);
			neigbor.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 0, true);
			if(isFeasible(neigbor)){
				neigbor.aperturesToIntensityMatrixPerStationwithoutIntensty();
				neigbor.ApertureMatrixStationToIntensityMatrixPlan();
				neigbor.IntensityMatrixPlan();
				GS.reset();
				GS.setNewModel(neigbor);
				neigbor.setIntensity(GS.newIntensity);
				neigbor.setOF(GS.objVal);
				//System.out.println(GS.objVal);
				count++;
				if(neigbor.singleObjectiveValue<sol.singleObjectiveValue) {
					neighborhood.add(neigbor);
					System.out.println(count);
					return neighborhood;
				}
			}
			
			TreatmentPlan neigbor1=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor1.updateSolDAO(sol);
			neigbor1.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 1, true);
			if(isFeasible(neigbor1)){
				neigbor1.aperturesToIntensityMatrixPerStationwithoutIntensty();
				neigbor1.ApertureMatrixStationToIntensityMatrixPlan();
				neigbor1.IntensityMatrixPlan();
				GS.reset();
				GS.setNewModel(neigbor1);
				neigbor1.setIntensity(GS.newIntensity);
				neigbor1.setOF(GS.objVal);
				count++;
				if(neigbor1.singleObjectiveValue<sol.singleObjectiveValue) {
					neighborhood.add(neigbor1);
					System.out.println(count);
					return neighborhood;
				}
			}
			TreatmentPlan neigbor2=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor2.updateSolDAO(sol);
			neigbor2.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 0, false);
			if(isFeasible(neigbor2)){
				neigbor2.aperturesToIntensityMatrixPerStationwithoutIntensty();
				neigbor2.ApertureMatrixStationToIntensityMatrixPlan();
				neigbor2.IntensityMatrixPlan();
				GS.reset();
				GS.setNewModel(neigbor2);
				neigbor2.setIntensity(GS.newIntensity);
				neigbor2.setOF(GS.objVal);
				count++;
				if(neigbor2.singleObjectiveValue<sol.singleObjectiveValue) {
					neighborhood.add(neigbor2);
					System.out.println(count);
					return neighborhood;
				}
			}
			TreatmentPlan neigbor3=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			neigbor3.updateSolDAO(sol);
			neigbor3.apertures.get(angleChoose).get(apertureChoose).changeRow(i, 1, false);
			//neigbor.changeAperture(aperture, angleChoose, apertureChoose);
			if(isFeasible(neigbor3)){
				neigbor3.aperturesToIntensityMatrixPerStationwithoutIntensty();
				neigbor3.ApertureMatrixStationToIntensityMatrixPlan();
				neigbor3.IntensityMatrixPlan();
				GS.reset();
				GS.setNewModel(neigbor3);
				neigbor3.setIntensity(GS.newIntensity);
				neigbor3.setOF(GS.objVal);
				count++;
				if(neigbor3.singleObjectiveValue<sol.singleObjectiveValue) {
					neighborhood.add(neigbor3);
					System.out.println(count);
					return neighborhood;
				}
			}
	
		}

		return neighborhood;
	}
	
	
	
	//_---------------------- Funciones intermedias para Algoritmos-------------------------//
	
	/**
	 * Esta funcion retorna el mejor vecino obtendio de un movimiento 
	 * @param neighborhood Lista de vecinos obtenidos a traves de un movimiento
	 * @param function que funncion se va a utilizar para saber que vecino es mejor
	 * @return Mejor vecino
	 */
	public static TreatmentPlan bestNeigbor(ArrayList<TreatmentPlan> neighborhood,int function,DDM M) {

		 
		 
		 
		TreatmentPlan best=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		best.singleObjectiveValue=9999999;
		for(TreatmentPlan actual:neighborhood) {
			actual.aperturesToIntensityMatrixPerStation();
			actual.intensityMatrixStationToIntensityMatrixPlan();
			evaluateSolution(actual,M,o,zmin,zmax,actual.weights);
			if(actual.singleObjectiveValue<best.singleObjectiveValue) {
				best.updateSolDAO(actual);
				best.aperturesToIntensityMatrixPerStation();
				best.intensityMatrixStationToIntensityMatrixPlan();
				evaluateSolution(best,M,o,zmin,zmax,best.weights);
			}
		}
		
		
		return best;
	}
	
	public static TreatmentPlan bestNeigborSolver(ArrayList<TreatmentPlan> neighborhood,int function,DDM M,Gurobi_Solver GS) {

		 
		 
		
		int count=0;
		TreatmentPlan best=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		best.singleObjectiveValue=9999999;
		for(TreatmentPlan actual:neighborhood) {
			actual.aperturesToIntensityMatrixPerStationwithoutIntensty();
			actual.ApertureMatrixStationToIntensityMatrixPlan();
			actual.IntensityMatrixPlan();
			GS.reset();
			GS.setNewModel(actual);
			actual.setIntensity(GS.newIntensity);
			actual.setOF(GS.objVal);
			if(actual.singleObjectiveValue<best.singleObjectiveValue) {
				best.updateSolDAO(actual);
				
			}
			count++;
		}
		System.out.println(count++);
		
		return best;
	}

	public static void updateFrontier(ArrayList<TreatmentPlan> neighborhood,ArrayList<TreatmentPlan> actualFront,DDM M) {
		boolean domino;
        //System.out.println("Entre a la funcion");
        ArrayList<TreatmentPlan> dominados=new ArrayList<TreatmentPlan>();//Se guardan en una lista para no generar problemas
        for(TreatmentPlan  neighbor: neighborhood ) {
        	//System.out.println("Neigbor");
        	//neighbor.datosFitness();
        	//System.out.println("");
        	
        	neighbor.aperturesToIntensityMatrixPerStation();
        	neighbor.intensityMatrixStationToIntensityMatrixPlan();
			evaluateSolution(neighbor,M,o,zmin,zmax,neighbor.weights);
            domino=false;
            for(TreatmentPlan nodo : actualFront ) {
            	nodo.aperturesToIntensityMatrixPerStation();
            	nodo.intensityMatrixStationToIntensityMatrixPlan();
    			evaluateSolution(nodo,M,o,zmin,zmax,nodo.weights);
                if(nodo.isDominatedEvaluetion(neighbor)) { 
                	//System.out.print("entre");
                	//nodo.print();
                	//neighbor.print();
                    domino=true;
                    dominados.add(nodo);
                    
                }

            }
            if(domino) {
            	//System.out.println("Si domino");
            	actualFront.add(neighbor);
                if(!dominados.isEmpty()) {
                	actualFront.removeAll(dominados);
                }
                dominados.clear();
            }
            
        }
        neighborhood.clear();
	}
	
	
	
	
	//---------------------------Funciones de Lectura----------------------------------------//
	
	public static void readInputFile(String dir) throws IOException{
	        
    	String sp="\\s+";
        //String dir = "./inputFile.txt";
        File f = new File(dir);
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line=fileIn.readLine();
        //First read numbero of Organs and angles
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                numOrgans = Integer.parseInt(auxReader[0]);
                numAngles = Integer.parseInt(auxReader[1]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        //Go to the next input line
        while(line != null){
            if (!line.contains("%")){
                break;
            }
            line=fileIn.readLine();
        }
 
        //Info Organs
        o = new Organs[numOrgans];
        while(line != null){
            if (!line.contains("%")){
                for (int y=0;y<numOrgans;y++){
                    String[] auxReader = line.split(sp);
                    o[y]=new Organs(
                        auxReader[0], 
                        Integer.parseInt(auxReader[1]), 
                        Integer.parseInt(auxReader[2]),
                        Integer.parseInt(auxReader[3]),
                        Integer.parseInt(auxReader[4]),
                        Integer.parseInt(auxReader[5]),
                        Integer.parseInt(auxReader[6]),
                        Integer.parseInt(auxReader[7]),
                        Integer.parseInt(auxReader[8]),
                        Integer.parseInt(auxReader[9]), 
                        Integer.parseInt(auxReader[10]),
                        Integer.parseInt(auxReader[11]), 
                        Integer.parseInt(auxReader[12]), 
                        Boolean.parseBoolean(auxReader[13]));
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        //get filepath
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                pathFile=auxReader[0];
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get option
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                option=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get selectionCriterion (1 = LS, 2 = nextDescent, 3= gradientBas)
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                selectionCriterion=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get max iterations LS
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                max_iter=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get max time LS
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                max_time=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get Vx Index. % of the volume receiving more than x% of the prescribed dose
        Vx = new int[numOrgans][];
        while(line != null){
            if (!line.contains("%")){
                while(!line.contains("%")){
                    String[] auxReader = line.split(sp);
                    Vx[Integer.parseInt(auxReader[0])]=new int[auxReader.length-1];
                    for(int i=1;i<auxReader.length;i++){
                        Vx[Integer.parseInt(auxReader[0])][i-1]=Integer.parseInt(auxReader[i]);
                    }
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        //get Dx Index. Minimum dose of the hottest x% of the volume
        Dx = new int[numOrgans][];
        while(line != null){
            if (!line.contains("%")){
                while(!line.contains("%")){
                    String[] auxReader = line.split(sp);
                    Dx[Integer.parseInt(auxReader[0])]=new int[auxReader.length-1];
                    for(int i=1;i<auxReader.length;i++){
                        Dx[Integer.parseInt(auxReader[0])][i-1]=Integer.parseInt(auxReader[i]);
                    }
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        
        //get initial BACs
        initialBACs = new int[numAngles];
        while(line != null){
            if (!line.contains("%")){
                while(!line.contains("%")){
                	String[] auxReader = line.split(sp);
                    for (int i=0; i<numAngles;i++){
                        initialBACs[i]= Integer.parseInt(auxReader[i]);
                    }
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        
        //Initial value aperture intensity
        while(line != null){
            if (!line.contains("%")){
            	String[] auxReader = line.split(sp);
            	init_intensity=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        //Maximum intensity for apertures
        while(line != null){
            if (!line.contains("%")){
            	String[] auxReader = line.split(sp);
            	max_intensity=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        //Máximo delta para la variación de intensidad en las aperturas
        while(line != null){
            if (!line.contains("%")){
            	String[] auxReader = line.split(sp);
            	max_delta=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
      //Step size for aperture intensity
        while(line != null){
            if (!line.contains("%")){
            	String[] auxReader = line.split(sp);
            	step_intensity=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        fileIn.close();
	}
	 
	//-------------------------Funciones de evaluacion-------------------------------------------------//
	public static double evaluateSolution(TreatmentPlan solution, DDM M, Organs[] o, int[] Zmin, int[] Zmax, double w[]){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		
		double[] solutionVector = solution.intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		//double[] solutionVector = new double[]{4,4,4,21,21,15,10,4,4,21,21,21,21,15,15,10,10,4,21,21,15,15,15,15,10,8,4,21,15,15,15,15,10,10,10,4,21,15,15,10,10,10,10,10,4,21,15,15,15,15,10,10,10,4,15,21,21,15,10,10,10,4,10,21,10,8,4,4,4,4,10,10,10,11,4,4,4,10,10,10,17,17,14,11,4,10,10,10,17,14,14,14,4,11,11,11,11,11,17,10,4,4,10,10,14,14,14,11,10,4,4,14,14,14,14,14,14,4,4,10,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,9,4,4,4,4,4,4,9,9,9,9,9,4,4,4,9,13,13,13,14,6,4,4,4,4,4,6,6,9,9,4,4,4,4,4,4,4,9,9,4,4,4,4,4,4,4,14,4,4,4,4,9,9,9,9,9,4,4,6,9,6,4,4,4,4,4,4,4,4,7,5,4,4,5,9,9,9,4,4,7,7,13,13,9,5,4,4,5,13,13,13,13,9,9,5,4,9,9,9,9,13,9,9,9,4,7,9,9,9,9,9,13,5,4,7,7,7,7,7,7,7,9,4,5,7,7,7,7,7,9,13,4,4,4,5,7,4,4,4,4,4,4,4,4,4,4,4,4,7,4,4,4,12,11,11,11,11,11,7,4,4,7,21,12,12,12,12,12,4,4,4,21,12,12,12,7,7,7,4,4,21,12,12,12,12,12,4,4,11,11,11,11,11,12,12,4,7,11,11,11,12};
		//double[] solutionVector = new double[]{11,20,18,18,20,20,18,18,4,18,18,18,18,18,18,20,18,4,11,18,18,18,18,18,20,18,18,11,18,18,18,18,18,20,18,18,11,18,18,20,18,18,18,18,18,18,20,20,18,18,18,18,18,18,20,20,20,20,18,18,18,4,20,20,18,11,11,4,11,11,10,16,16,10,10,4,10,10,10,10,10,10,10,13,16,16,16,16,16,13,13,13,16,16,16,16,16,16,13,13,4,16,16,16,13,13,13,13,13,4,13,13,16,16,13,13,13,10,4,10,10,13,13,13,10,10,4,4,4,13,4,4,10,4,4,4,4,11,4,4,4,5,5,5,11,8,8,8,4,5,5,5,11,8,8,8,4,4,5,5,5,11,8,8,8,4,4,5,5,5,11,8,8,8,4,4,5,5,5,11,8,8,8,4,4,5,5,8,11,8,8,5,4,8,8,8,8,4,4,8,8,4,4,4,5,4,4,4,5,5,5,4,4,4,5,5,5,5,5,9,4,4,14,5,5,5,5,5,4,4,4,14,5,5,5,5,5,4,4,4,14,5,5,5,4,4,4,4,4,14,5,5,5,5,5,4,4,4,14,9,9,9,9,9,4,4,4,5,9,5,4,4,4,4,4,4,4,4,4,4,9,8,0,0,0,0,0,9,9,9,8,6,4,4,4,4,9,9,9,9,6,4,4,4,4,9,9,9,9,6,4,4,4,4,9,9,9,9,6,4,4,4,9,8,8,8,8,4,4,4,4,4,8,4,4,4};
		double F = 0.0, pen,score;
		
		//System.out.println("Largo solución :"+solutionVector.length);
		//Recorremos los organos
		for(int i=0; i<o.length; i++){
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			//Recorremos claves de voxel por organo para su evaluación
			while(keys.hasMoreElements()){
	            key = keys.nextElement();
	            beams = aux_index.get(key);
	            intensityVoxel2Beams = 0.0;
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation;
	         
	            }
	            
	            /*
	            if(intensityVoxel2Beams < Zmin[i]){
            		pen += w[i] * (Math.pow(Zmin[i]-intensityVoxel2Beams, 2) );
            	}
   				if(intensityVoxel2Beams > Zmax[i] ){
   					pen += w[i] * (Math.pow(intensityVoxel2Beams-Zmax[i], 2) );
   				}
   				*/
            	if(i == 0){
            		pen += w[i] * Math.pow(Math.max((Zmin[i] - intensityVoxel2Beams),0),2);
            	}else{
   					pen += w[i] * Math.pow(Math.max((intensityVoxel2Beams - Zmax[i]),0),2);
   				}
            	
	        }
			//System.out.println(i+" "+Zmax[i]+" ");
			score=pen/aux_index.size();
			solution.scores[i]=score;
			F+=pen/aux_index.size();
		}
		solution.singleObjectiveValue=F;
		//System.out.println("Solucion: "+F);
		return F;
	}
	
	/**
	 * Esta funcion retorna verdadero si es que la solucion es factible, en caso contrario devuelve falso
	 * @param sol
	 * @return
	 */
	public static boolean isFeasible(TreatmentPlan sol){
		Vector<List<Aperture>> stations = sol.apertures;
		boolean feasible=true;
		//recorrido de angulos
		for(int i=0; i<stations.size(); i++){
			
			List<Aperture> list_apertures = stations.get(i);
			int[][] limit_shape = initial_aperture_shape.get(i);
			//recorrido de aperturas por angulo
			for(int j=0; j<list_apertures.size(); j++){
				//int intensity = list_apertures.get(j).intensity;
				int[][] matrix = list_apertures.get(j).aperture;
				
				for(int x=0;x<matrix.length; x++){
					if(limit_shape[x][0]>matrix[x][0] || limit_shape[x][1]<matrix[x][1]) {
						//System.out.println(limit_shape[x][0]+ " " + matrix[x][0]+ " "+limit_shape[x][1]+ " " + matrix[x][1]);
						feasible=false;
						return feasible; 
					}
					if(matrix[x][0]>matrix[x][1]) {
						//System.out.println(limit_shape[x][0]+ " " + matrix[x][0]+ " "+limit_shape[x][1]+ " " + matrix[x][1]);
						feasible=false;
						return feasible; 
					}

					
				}	
			}
				
		}
		return feasible;
	}
	
	public static void setDaoExample(TreatmentPlan sol) {
		
	int [][][]aper= {{{0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,1,0,1,1,1,0,1,0,1,1,0,1,0,1,0,0,0,1,1,0,0,1,0,1,0,1,0,0,1,0,1,1,1,0,1,1,0,0,1,0,0,0,1,0,1,0,0,0
	},{1,0,0,0,0,1,0,1,1,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,1,1,1,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,0,1
	},{0,0,1,1,0,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,0,0,1,0,1,0,0,1,1,0,0,1,0,0,1,0,0,1,0,1,1,1,0,0,1,1,0,1,1,1,0,1,1,0,0,0,1
	},{1,1,0,1,0,0,0,1,0,1,0,1,1,1,0,0,0,0,1,0,1,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,0,1,1,0,0,1,0,1,1,0,0,0,0
	},{1,1,0,1,0,1,1,1,0,1,1,0,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,0,0,1,1,1,0,0,0,1,0,1,1,0,1,0,0,0,1,1,1,1,1,0,0,1,1,1,0,1,1,1,1,1,0,1,0,1,0,0,1,1,0
	}},{{0,0,1,0,0,0,0,1,1,1,1,0,1,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,1,1,1,1,1,1,0,1,0,0,0,0,1,1,1,0,1,0,0,1,0,1},{1,1,1,1,1,1,0,0,1,1,1,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,0,1,0,0,0,0,1,0,0,0,1,0,0,1,1,1,1,0,1,0,1,1,0,0,0,1},{0,1,0,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,0,1,0,1,0,0,1,1,1,0,1,1,0,0,0,1,0,1,1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1},{1,1,1,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,1,0,0,1,0,0,0,0,0,1,0,1,1,1,0,0,1,1,0,1,0,1,0,0,1,0,0,1,0,0,1,0,1,1,0,0,0,1,1,0,0,1,1,1,0,1,1},{0,1,0,1,0,1,1,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,0,1,1,1,1,0,1,0,1,1,1,1,0,1,1,0,0,1,1,1,0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,1,0,1,1}
	},{{0,0,0,1,0,1,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,1,0,1,1,0,1,1,1,0,0,1,1,1,1,0,0,0,1,1,0,0,1},{0,1,1,0,0,1,0,1,0,0,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,0,1,0,1,1,1,0,0,1,0,0,0,1,0,0,1,1,1,0,0,1,1,1,1,0,0,1,0,0,1,1,0,1,0,0,0,0,1,0,0,1,1},{0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,1,0,0,1,1,0,1,1,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,1,0,1,1,0,1,0,0,0,1,1,1,1,0,1,0,1},{1,0,1,0,1,0,1,1,0,1,0,1,1,1,0,0,1,1,0,0,0,1,1,1,1,0,1,1,0,1,0,1,0,0,1,0,0,1,1,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,1,0,0,1,1,1,1,0,0,1,0,1,0},{1,1,0,1,1,1,1,1,1,0,1,1,0,1,1,0,0,0,1,1,0,0,0,1,0,1,1,1,0,1,0,0,1,0,1,0,1,1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,1}},{{1,1,1,0,0,1,0,0,0,1,0,1,0,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,1,1,0,0,1,1,0,0,0,1,0,0,0},{1,0,0,1,1,1,0,0,1,0,1,0,0,1,0,1,0,1,1,0,0,0,0,0,1,1,0,0,0,1,1,1,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,1,1,1,0,1,0,1,0,0,1,1,1,1,1,0,1,1,1,0,0,1,1,0},{1,0,1,1,0,1,1,0,1,0,1,1,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,0,1,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,1,0,1,0,0,1,1},{0,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0,1,0,0,1,0,0,1,1,0,1,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,1,1,1,0,0,1,1,0,1,0,1,1,1,1,1,0,1,0,1,0,1,1,0},{0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,1,1,1,0,1,1,1,0,0,1,1,0,1,0,0,1,0,1,0,1,0,1,0,0,1,1,0,1,1,0,1,0,0,1,0,1,0,1,0,0,0,1,0,1,0,0,1,1,1,0,1}
	},{{0,1,0,1,0,1,0,0,0,1,0,1,0,0,1,1,0,1,0,1,0,1,0,0,1,1,0,1,0,1,1,0,0,0,1,0,1,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,1},{1,1,0,1,1,0,0,0,1,0,0,1,1,1,0,1,1,0,1,0,0,1,0,0,1,0,1,1,0,1,0,1,1,0,0,1,0,0,1,1,0,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,1,0,0,0,1},{0,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,0,1,1,0,0,1,0,1,0,0,0,0,0,1,0,1,1,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,1,1,1,1,0,1,1,0},{0,0,1,0,1,0,0,0,1,1,0,1,0,1,0,0,1,1,0,0,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,0,0,1,1,0,1,0,0,0,0,1,1,0,1,1,1,0,0,0,1,1,0,0,0,0,1},{0,0,1,0,1,1,1,0,1,0,0,1,0,1,1,0,0,0,0,1,1,0,0,1,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,0,1,0}
	}}; 	
		
	
	for(int i=0;i<sol.aperturesBmlts.size();i++) {
		double[][] intensities;
		intensities=sol.aperturesBmlts.get(i);
			for(int j=0;j<intensities.length;j++) {		
				for(int k=0;k<intensities[j].length;k++) {
					intensities[j][k]=aper[i][j][k];
				
				}
			}
		}
	}
	public static void setBeams(TreatmentPlan solution, String dir) throws IOException{
		totalbmlt=0;
		for(int i=0; i<numAngles; i++){
			solution.selAngles[i] = new Beam(beamletSet[initialBACs[i]][1],initialBACs[i], pathFile);
			totalbmlt=totalbmlt+solution.selAngles[i].beamlets;
		}
	}
	
	public static void setInitialApertureShape(TreatmentPlan solution){
		int sizeX, sizeY, x, y;
		apertureShapes=new Vector<int[][]>();
		initial_aperture_shape=new Vector<int[][]>();
		int beamlets=0;
		boolean leftLimit, rightLimit;
		
		for(int i=0;i<numAngles;i++){
			sizeX  = solution.selAngles[i].maxAbsX;
			sizeY = solution.selAngles[i].maxAbsY;
			int [][] apertureShape = new int[sizeX][sizeY];
			for(x=0; x<sizeX; x++){
				for(y=0; y<sizeY; y++){
					apertureShape[x][y] = -1;
				}
			}
			
			double[][] beamletsCoord = solution.selAngles[i].beamletsCoord;
			for(int j=0;j<solution.selAngles[i].beamlets;j++){
				int xcoord=(int)beamletsCoord[j][1]-1;
				int ycoord=(int)beamletsCoord[j][2]-1;
				apertureShape[xcoord][ycoord] = 0;
				beamlets=beamlets+1;	
			}
			apertureShapes.add(apertureShape);
			
			int[][] apertureVectorShape = new int[sizeX][2];
			for(x=0; x<sizeX; x++){ 
				leftLimit = false;
				rightLimit = false;
				for(y=0; y<sizeY; y++){
					if(!leftLimit && apertureShape[x][y] == 0){
						apertureVectorShape[x][0] = y;
						leftLimit = true;
					}
					if(leftLimit && !rightLimit && apertureShape[x][y] == 0){
						apertureVectorShape[x][1] = y;
						leftLimit = true;
					}
				}
			}
			
			initial_aperture_shape.add(apertureVectorShape);
			solution.beamlets=beamlets;
		}
	}
	
	//----------------------------Funciones para la generacion inicial de soluciones
	
	public static void generateFirstSolution(TreatmentPlan sol){
		int [][] aper_shape, new_aper;
		int range;
		Random r = new Random();
		int random_intensity;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			for(int z=0;z<num_apertures;z++){
				new_aper = new int[aper_shape.length][2];
				random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
				
				//Solución inicial con limites random
				for(int j=0; j<aper_shape.length; j++){
					range=aper_shape[j][1]-aper_shape[j][0];
					new_aper[j][0] = ((int) r.nextInt(range+1))+ aper_shape[j][0];
					range=aper_shape[j][1]-new_aper[j][0];
					new_aper[j][1] = (int) r.nextInt(range+1) + new_aper[j][0];
				}
				
				Aperture aper = new Aperture(new_aper, random_intensity);
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	
	
	public static void generateFirstSolutionWithFixedIntensity(TreatmentPlan sol){
		int [][] aper_shape, new_aper;
		int range;
		Random r = new Random();
		List<Aperture> list_aperture;
		sol.apertures.clear();
		
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			for(int z=0;z<num_apertures;z++){
				new_aper = new int[aper_shape.length][2];
				
				//Solución inicial con limites random
				for(int j=0; j<aper_shape.length; j++){
					range=aper_shape[j][1]-aper_shape[j][0];
					new_aper[j][0] = ((int) r.nextInt(range+1))+ aper_shape[j][0];
					range=aper_shape[j][1]-new_aper[j][0];
					new_aper[j][1] = (int) r.nextInt(range+1) + new_aper[j][0];
				}
				Aperture aper = new Aperture(new_aper, init_intensity);
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	public static void generateFirstSolutionOpen(TreatmentPlan sol){
		int [][] aper_shape, new_aper;

		Random r = new Random();
		int random_intensity;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			for(int z=0;z<num_apertures;z++){
				new_aper = new int[aper_shape.length][2];
				random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
				
				
				//Solución inicial completamente abierta
				for(int j=0; j<aper_shape.length; j++){
					new_aper[j][0] = aper_shape[j][0];
					
					new_aper[j][1] = aper_shape[j][1];
				}
				
				Aperture aper = new Aperture(new_aper, random_intensity);
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	
	public static void generateFirstSolutionOpenWithFixedIntensity(TreatmentPlan sol){
		int [][] aper_shape, new_aper;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			for(int z=0;z<num_apertures;z++){
				new_aper = new int[aper_shape.length][2];

				//Solución inicial completamente abierta
				for(int j=0; j<aper_shape.length; j++){
					new_aper[j][0] = aper_shape[j][0];
					new_aper[j][1] = aper_shape[j][1];
				}
				
				Aperture aper = new Aperture(new_aper, init_intensity);
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	
	
	
	//--------------------------- Funciones para imprimir la solucion
	public static void printFront(ArrayList<TreatmentPlan> Front) {
		for(TreatmentPlan nodo:Front) {
			nodo.print();
		}
	}

	public static void printFirstSolution(TreatmentPlan sol){
		Vector<List<Aperture>> stations = sol.apertures;
		
		//recorrido de angulos
		for(int i=0; i<stations.size(); i++){
			System.out.println("Angulo: "+initialBACs[i]);
			List<Aperture> list_apertures = stations.get(i);
			//recorrido de aperturas por angulo
			for(int j=0; j<list_apertures.size(); j++){
				int intensity = (int) list_apertures.get(j).intensity;
				int[][] matrix = list_apertures.get(j).aperture;
				System.out.println("Apertura N°: "+j);
				System.out.println("Intensidad: "+intensity+"\n");
				
				for(int x=0;x<matrix.length; x++){
					System.out.println(matrix[x][0]+", "+matrix[x][1]);
				}
				
			}
			System.out.println("\n");
		}
		
	}
	public static void printIntensitiesMatrix(TreatmentPlan sol) {
		int sizeX,sizeY,x,y;
		for(int i=0;i<sol.intesities.size();i++) {
			sizeX  = sol.selAngles[i].maxAbsX;
			sizeY = sol.selAngles[i].maxAbsY;
			double[][] intensity;
			intensity=sol.intesities.get(i);
			for(x=0;x<sizeX;x++) {
				for(y=0;y<sizeY;y++) {
					System.out.print(intensity[x][y]+" ");
				}
				System.out.println();
			}
			System.out.println();
			
		}
	}
	
	public static void printapertureBeamlet(TreatmentPlan sol) {
		int sizeX,sizeY,x,y;
		for(int i=0;i<sol.aperturesBmlts.size();i++) {
			System.out.println("Angulosb"+i);
			for(x=0;x<sol.aperturesBmlts.get(i).length;x++) {
				for(y=0;y<sol.aperturesBmlts.get(i)[x].length;y++) {
					System.out.print(sol.aperturesBmlts.get(i)[x][y]+" ");
				}
				System.out.println();
			}
			System.out.println();
			
		}
	}
	
	public static void printShapeSolution(TreatmentPlan sol){

		//recorrido de angulos
		for(int i=0; i<initial_aperture_shape.size(); i++){

			System.out.println("Angulo: "+initialBACs[i]);
			int[][] list_shape = initial_aperture_shape.get(i);
			//recorrido de aperturas por angulo
	
				for(int x=0;x<list_shape.length; x++){
					System.out.println(list_shape[x][0]+", "+list_shape[x][1]);
				}
				
			
			System.out.println("\n");
		}
		
	}

}
