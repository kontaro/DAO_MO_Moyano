package IMRT_Base;

/* LogFunction
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Vector;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author gcab623
 * 
 * @Modificaciones Denisse Katherine Maicholl Mauricio
 */
public class TreatmentPlan {

    public int beams;
    public int beamlets;
    public double slope;
    public double[] intensity;
    public Beam[] selAngles;
    public double[] weights;
    public double[] gEUD;
    public double singleObjectiveValue;
    public boolean visited;
    public double Vx[][];
    public double Dx[][];
    public double scores[];//katty
    //public double w[];//katty se ocupa o[i].weight
    public double globalScore;//katty
    //public double Dp[];//katty
    public double Dxx[][];//katty dosis por cada organo,PTV D95 ;OAR Dmean Dmax D50 D75 D100
    public double Vxx[][];//katty
    public double allScores[][];//katty scores por cada organo
    public double quadFunction[]; //Maiik FO por cada organo
    public double epsilonDominance[]; // % of the objective value that is used to relax/stress the dominance condition. It might be positive or negative
    
    public Vector<int[][]> aperturesLimit;   // # DAO : vector con beam por angulo
    public Vector<double[][]> intesities;       // # DAO : vector de intesidades
    public Vector<List<Aperture>> apertures; // # DAO : vector para cada estaci髇 que dentro guarda ua lista de las N aperturas por cada una de ellas
    public int init_intensity;               // # DAO : Intesidad inicial para las aperturas
    public int max_intesity;                 // # DAO : Intesidad m醝ma para las aperturas
    public int max_delta;                    // # DAO : Delta m醲imo para la variaci髇 de intesidad por apertura
    public int max_iter;                     // # DAO : Cantidad m醲ima de iteraciones del algoritmo
    public int max_time;                     // # DAO : tiempo m醲imo de ejecuci髇 del algoritmo
    public int seed;                         // # DAO : semilla para partir la busqueda local
    public int step_intensity;               // # DAO : Step size for aperture intensity (2)
    public Vector<double[][]>aperturesBmlts;
    
    public TreatmentPlan(int b, int bl, int a, int o) {
        this.beams = b;
        this.beamlets = bl;
        this.intensity = new double[bl];
        this.selAngles = new Beam[a];
        this.gEUD = new double[o + 1]; // # of Organs + PTV UB
        this.weights = new double[o]; // # of Organs + PTV UB
        this.Vx = new double[o][]; // # of Organs
        this.Dx = new double[o][]; // # of Organs
        this.singleObjectiveValue = 0;
        this.visited = false;
        this.slope = 0;
        this.epsilonDominance = new double[o];//katty
        this.scores = new double[o];//katty
        //this.w = new double[o];//katty
        this.globalScore = 0;//katty
        //this.Dp = new double[o];//katty
        this.Dxx = new double[o][];//katty
        this.Vxx = new double[o][];//katty
        this.allScores = new double [o][];//katty
        this.quadFunction = new double[o];//Maicholl
    }
    
    /************************************************************
	 ************** Constructor sobrecargado para DAO ************
	 *************************************************************/
	 public TreatmentPlan(int init, int max_i, int delta, int iter, int time, int seed, int step, int angles,int numOrgans) {
	 	this.beams = angles;
	 	this.selAngles = new Beam[angles];
	 	this.intesities = new Vector<double[][]>();
	 	this.apertures = new Vector<List<Aperture>>();
	 	this.aperturesBmlts = new Vector<double[][]>();
	 	this.aperturesLimit=new Vector<int[][]>();
	 	this.init_intensity = init;
	 	this.max_intesity = max_i;
	 	this.max_delta = delta;
	 	this.max_iter = iter;
	 	this.max_time = time;
	 	this.seed = seed;
	 	this.step_intensity = step;
	 	this.singleObjectiveValue = 0;
	 	this.scores=new double[numOrgans];
	 	
	 	
	}
	public void copyMatrix(Vector<int[][]> Original,Vector<int[][]> copy) {
		copy.clear();

		for(int[][] auxAperture:Original) {
			int [][]aux= new int[auxAperture.length][auxAperture[0].length];
			for(int i=0;i<auxAperture.length;i++) {
				for(int j=0;j<auxAperture[0].length;j++) {
					aux[i][j]= auxAperture[i][j];
				}
			}
			copy.add(aux);
		}
		
	
	}
	public double[][] copyMatrix(int[][] Original) {

		  double [][]copy= new double[Original.length][Original[0].length];
			for(int i=0;i<Original.length;i++) {
				for(int j=0;j<Original[0].length;j++) {
					copy[i][j]= Original[i][j];
				}
			}
			return copy;
		
	
	}
	public void aperturesToIntensityMatrixPerStation(){
    	Vector<List<Aperture>> vectorStations = apertures;
    	//Vector<int[][]> vectorIntensities = tp.intesities;
    	intesities.clear();
    	//recorrido de estaciones
    	for(int i=0; i<vectorStations.size(); i++){
    		List<Aperture> listApertures = vectorStations.elementAt(i);  
    		double[][] intensitiesMatrix = copyMatrix(aperturesLimit.get(i));
    		
    		//recorrido de lista de aperturas por estaci髇
            for(int j=0; j<listApertures.size(); j++ ){
            	int[][] aperture = listApertures.get(j).aperture;
            	double intensity = listApertures.get(j).intensity;
            	
            	//recorrido de fila de matriz de intensidades
            	for(int row=0; row<intensitiesMatrix.length; row++){
            		int leftLimit = aperture[row][0];
            		int rightLimit = aperture[row][1];
            		
            		//recorrido de columnas de matrix asinando intensidades seg鷑 limites
            		for(int column=leftLimit; column<=rightLimit; column++){
                		
            			intensitiesMatrix[row][column] =intensitiesMatrix[row][column] + intensity;
            			
                	}
            	}
            }
            intesities.add(intensitiesMatrix);
        }
    	
    }
	public void aperturesToIntensityMatrixPerStationwithoutIntensty(){
    	Vector<List<Aperture>> vectorStations = apertures;
    	//Vector<int[][]> vectorIntensities = tp.intesities;
    	intesities.clear();
    	//recorrido de estaciones
    	for(int i=0; i<vectorStations.size(); i++){
    		List<Aperture> listApertures = vectorStations.elementAt(i);  
    		
    		
    		//recorrido de lista de aperturas por estaci髇
            for(int j=0; j<listApertures.size(); j++ ){
            	int[][] aperture = listApertures.get(j).aperture;
            	double[][] intensitiesMatrix = copyMatrix(aperturesLimit.get(i));
            	
            	//recorrido de fila de matriz de intensidades
            	for(int row=0; row<intensitiesMatrix.length; row++){
            		int leftLimit = aperture[row][0];
            		int rightLimit = aperture[row][1];
            		
            		//recorrido de columnas de matrix asinando intensidades seg鷑 limites
            		for(int column=leftLimit; column<=rightLimit; column++){
                		
            			intensitiesMatrix[row][column] =intensitiesMatrix[row][column] + 1;
            			
                	}
            		
            	}
            	intesities.add(intensitiesMatrix);
            }
            
        }
    	
    }
	public void changeAperture(int[][] newAperture,int angle, int positionAperture) {
		apertures.get(angle).get(positionAperture).aperture=newAperture;
		//apertures.get(angle).add(newAperture);
		
	}
    public int[] intensityMatrixStationToIntensityMatrixPlan(Vector<int[][]> stationIntensities){
    	int numColumns = (stationIntensities.size())*(stationIntensities.elementAt(0).length)*(stationIntensities.elementAt(0)[0].length);
    	int [] solutionVector = new int[numColumns];
    	int columns = 0;
    	
    	//recorrido de estaciones
    	for(int i=0; i<stationIntensities.size(); i++){
    		int[][] intensityAperture = stationIntensities.elementAt(i);
    		
    		for(int row=0; row < intensityAperture.length; row++){
    			for(int column=0; column < intensityAperture[row].length; column++){
    				solutionVector[columns] = intensityAperture[row][column];
        		}
    		}
    	}
    	return solutionVector;
    }
    public int[] bealetsForAngle() {
    	int[] aux=new int[selAngles.length];
    	for(int j=0;j<aux.length;j++){
			aux[j]=selAngles[j].beamlets;
		}
    	
    	return  aux;
    }
    public void intensityMatrixStationToIntensityMatrixPlan(){
    	
    	intensity=new double[beamlets] ;
    	int index=0;
    	for(int i=0;i<intesities.size();i++){
    		double[][] intensities;
    		intensities=intesities.get(i);
			double[][] beamletsCoord = selAngles[i].beamletsCoord;
			for(int j=0;j<selAngles[i].beamlets;j++){
				int xcoord=(int)beamletsCoord[j][1]-1;
				int ycoord=(int)beamletsCoord[j][2]-1;
				//System.out.print(xcoord+" ");
				//System.out.println(ycoord);
				//System.out.println(beamlets);
				//System.out.println(intensities.length);
				intensity[index]=intensities[xcoord][ycoord];
				index=index+1;
			}
    	}
    	//recorrido de estaciones

    	
    }
    public void IntensityMatrixPlan(){
    	
    	
    	
    	
    	intensity=new double[beamlets] ;
    	for(int t=0;t<beamlets;t++) {intensity[t]=0;}
    	int index=0;
    	for(int a=0;a<aperturesBmlts.size();a++) {
	    	double[][] intensities;
	    	
	    	for(int i=0;i<aperturesBmlts.get(a).length;i++){
	    		intensities=aperturesBmlts.get(a);
				
				for(int j=0;j<selAngles[a].beamlets;j++){
					
					//System.out.print(xcoord+" ");
					//System.out.println(ycoord);
					//System.out.println(beamlets);
					//System.out.println(intensities.length);
					intensity[(j+index)]=intensity[(j+index)]+intensities[i][j]*apertures.get(a).get(i).intensity;
					
				}
				
	    	}
	    	index=index+selAngles[a].beamlets;
	    	//System.out.println(index);
    	}
    	

    	
    }
    public void ApertureMatrixStationToIntensityMatrixPlan(){
    	
    	
    	int index=0;
    	for(int x=0;x<beams;x++) {
    		double[][] newAperture = new double[apertures.get(x).size()][selAngles[x].beamlets];
    		
	    	for(int i=0;i<apertures.get(x).size();i++){
	    		
	    		double[][] intensities;
	    		intensities=intesities.get(index);
	    		
				double[][] beamletsCoord = selAngles[x].beamletsCoord;
				for(int j=0;j<selAngles[x].beamlets;j++){
					int xcoord=(int)beamletsCoord[j][1]-1;
					int ycoord=(int)beamletsCoord[j][2]-1;
					newAperture[i][j]=intensities[xcoord][ycoord];
					
				}
				index=index+1;
	    	}
	    	aperturesBmlts.add(newAperture);
	    	
    	}
    	//recorrido de estaciones

    	
    }
    public void setWeights(double weight[]) {
    	this.weights=weight.clone();
    }
    
    
    public void updateSol(TreatmentPlan s) {
        this.beamlets = s.beamlets;
        this.beams = s.beams;
        this.singleObjectiveValue = s.singleObjectiveValue;
        this.intensity = new double[this.beamlets];
        this.selAngles = new Beam[this.beams];
        this.visited = s.visited;
        this.gEUD = new double[s.gEUD.length];
        this.epsilonDominance=s.epsilonDominance;
        this.globalScore = s.globalScore;//katty
        System.arraycopy(s.epsilonDominance, 0, this.epsilonDominance, 0, s.epsilonDominance.length);//katty
        System.arraycopy(s.intensity, 0, this.intensity, 0,
                Math.min(s.intensity.length, this.intensity.length));
        System.arraycopy(s.selAngles, 0, this.selAngles, 0, this.beams);
        System.arraycopy(s.gEUD, 0, this.gEUD, 0, s.gEUD.length);
        System.arraycopy(s.weights, 0, this.weights, 0, s.weights.length);
        //System.arraycopy(s.w, 0, this.w, 0, s.w.length);//katty
        System.arraycopy(s.scores, 0, this.scores, 0, s.scores.length);//katty
        System.arraycopy(s.quadFunction, 0, this.quadFunction, 0, s.quadFunction.length);//Maicholl
        //System.arraycopy(s.Dp, 0, this.Dp, 0, s.Dp.length);//katty

        for (int i = 0; i < s.Dxx.length; i++) {//katty
            if (null != s.Dxx[i]) {
                this.Dxx[i] = new double[s.Dxx[i].length];
                System.arraycopy(s.Dxx[i], 0, this.Dxx[i], 0,
                        Math.min(s.Dxx[i].length, this.Dxx[i].length));
            }
        }
        for (int i = 0; i < s.allScores.length; i++) {//katty
            if (null != s.allScores[i]) {
                this.allScores[i] = new double[s.allScores[i].length];
                System.arraycopy(s.allScores[i], 0, this.allScores[i], 0,
                        Math.min(s.allScores[i].length, this.allScores[i].length));
            }
        }
        for (int i = 0; i < s.Dx.length; i++) {
            if (null != s.Dx[i]) {
                this.Dx[i] = new double[s.Dx[i].length];
                System.arraycopy(s.Dx[i], 0, this.Dx[i], 0,
                        Math.min(s.Dx[i].length, this.Dx[i].length));
            }
        }
        for (int i = 0; i < s.Vx.length; i++) {
            if (null != s.Vx[i]) {
                this.Vx[i] = new double[s.Vx[i].length];
                System.arraycopy(s.Vx[i], 0, this.Vx[i], 0,
                        Math.min(s.Vx[i].length, this.Vx[i].length));
            }
        }
        this.slope = s.slope;
    }
    
    public void updateSolDAO(TreatmentPlan s) {
  	
    	
    	
        this.beamlets = s.beamlets;
        this.beams = s.beams;
        this.selAngles = new Beam[this.beams];
        this.singleObjectiveValue = s.singleObjectiveValue;
        this.init_intensity = s.init_intensity;
	 	this.max_intesity = s.max_intesity;
	 	this.max_delta = s.max_delta;
	 	this.max_iter = s.max_iter;
	 	this.max_time = s.max_time;
	 	this.seed = s.seed;
	 	this.step_intensity = s.step_intensity;
	 	this.weights=s.weights.clone();
	 	
	 	
        this.intensity = new double[this.beamlets];
	 	this.intesities = new Vector<double[][]>();
	 	
	 	this.apertures = new Vector<List<Aperture>>();
	 	copyMatrix( s.aperturesLimit, this.aperturesLimit);
	 	
	 	
	 	System.arraycopy(s.selAngles, 0, this.selAngles, 0, this.beams);
	 	//System.arraycopy(s.intensity,0,this.intensity,0,this.intensity.length);
	 	
	 	for (int i = 0; i < s.apertures.size(); i++) {
	 		List<Aperture> nuevo =  new ArrayList<Aperture>();
	 		for(int j=0;j<s.apertures.get(i).size();j++) {
	 			Aperture nuevaApertura=new Aperture(s.apertures.get(i).get(j));
	 			nuevo.add(nuevaApertura);
	 		}
	 		apertures.add(nuevo);
	 	}
	 	System.arraycopy(s.scores, 0, this.scores, 0, s.scores.length);
        

    }
    
    public void setIntensity(double[][] newIntensity) {
		for(int i=0; i<apertures.size(); i++){
			
			List<Aperture> list_apertures = apertures.get(i);
			//recorrido de aperturas por angulo
			for(int j=0; j<list_apertures.size(); j++){
				//list_apertures.get(j).intensity=Math.round(newIntensity[i][j]);
				list_apertures.get(j).intensity=newIntensity[i][j];
				
			}
			
		}
    }
    public void setOF(double newObjectiveValue) {
		this.singleObjectiveValue=newObjectiveValue;
    }
    

    public void generateReferencePoint(Organs[] o, int opt, String jobThreadID) 
            throws IOException {

        switch (opt) {
            case 1:
                updateWeights(o);
                this.solveLogFunction(o, jobThreadID);
                break;
            case 2:
                this.solveLexiRectum(o, jobThreadID);
                break;
            case 4:
                //this.solveUnconstrainedWeightedSum_gEUD(o, jobThreadID);
                break;
            case 5:
                this.solveConvexLogFunction(o, jobThreadID);
                break;
            case 6:
                this.solveWeightedSum(o, jobThreadID);
                break;
            case 7:
                this.solveAdaptiveWeightedSum(o, jobThreadID);
                break;
            case 8:
                //this.solveULogFunction(o, jobThreadID);
                break;
            case 9:
                this.solveCuadraticFunction(o,"ipopt", -1, jobThreadID, 0); //Maicholl
                break;
        }
        

    }
    
    /**
     *
     * @param o
     * @param solver
     * @param opt
     * @param jobThreadID
     * @param scoreConfigFile
     * @throws IOException
     */
    public void generateReferencePoint(Organs[] o, String solver, int opt, String jobThreadID, String scoreConfigFile) throws IOException {

        switch (opt) {
            case 1:
                updateWeights(o);
                this.solveLogFunction(o, solver, jobThreadID, scoreConfigFile);
                break;
            case 2:
                //this.solveLexiRectum(o,solver, jobThreadID, scoreConfigFile);
                break;
            case 4:
                updateWeights(o);
                this.solveUnconstrainedWeightedSum_gEUD(o,solver, jobThreadID, scoreConfigFile);
                break;
            case 5:
                //this.solveConvexLogFunction(o,solver, jobThreadID, scoreConfigFile);
                break;
            case 6:
                this.solveWeightedSum(o,solver, jobThreadID, scoreConfigFile);
                break;
            case 7:
                this.solveAdaptiveWeightedSum(o,solver, jobThreadID, scoreConfigFile);
                break;
            case 8:
                updateWeights(o);
                this.solveULogFunction(o, solver, jobThreadID, scoreConfigFile);
                break;
            case 9:
                updateWeights(o);
                this.solveCuadraticFunction(o,solver,20, jobThreadID,scoreConfigFile); //Maicholl
                break;
        }

    }
    
    public void generateReferencePoint(Organs[] o, String solver, int opt, String jobThreadID) throws IOException {

        switch (opt) {
            case 1:
                updateWeights(o);
                this.solveLogFunction(o, solver, jobThreadID);
                break;
            case 2:
                //this.solveLexiRectum(o,solver, jobThreadID, scoreConfigFile);
                break;
            case 4:
                updateWeights(o);
                this.solveUnconstrainedWeightedSum_gEUD(o,solver, jobThreadID);
                break;
            case 5:
                //this.solveConvexLogFunction(o,solver, jobThreadID, scoreConfigFile);
                break;
            case 6:
                this.solveWeightedSum(o,solver, jobThreadID);
                break;
            case 7:
                this.solveAdaptiveWeightedSum(o,solver, jobThreadID);
                break;
            case 8:
                updateWeights(o);
                this.solveULogFunction(o, solver, jobThreadID);
                break;
            case 9:
                updateWeights(o);
                this.solveCuadraticFunction(o,solver, 20, jobThreadID, 0); //Maicholl
                break;
        }

    }
    
    
    public void generateReferencePoint(Organs[] o, int weighted, String solver, double maxIntensity, int opt,  String jobThreadID) throws IOException {

        switch (opt) {
            case 1:
                updateWeights(o);
                this.solveLogFunction(o,solver, jobThreadID);
                this.calculateDoseAll(o,jobThreadID);//katy
                break;
            case 2:
                 this.solveLexiRectum(o,solver, jobThreadID);
                break;
            case 4:
                //this.solveUnconstrainedWeightedSum_gEUD(o,solver, jobThreadID);
                break;
            case 5:
                this.solveConvexLogFunction(o,solver, jobThreadID);
                break;
            case 6:
                this.solveWeightedSum(o,solver, jobThreadID);
                break;
            case 7:
                this.solveAdaptiveWeightedSum(o,solver, jobThreadID);
                break;
            case 8:
                updateWeights(o);
                this.solveULogFunction(o, solver, jobThreadID);
                break;
            case 9:
                this.solveCuadraticFunction(o,solver, maxIntensity, jobThreadID, weighted); //Maicholl
                break;
        }

    }
    /* //No debe existir ya que ya existen funciones del tipo organ, string,
    //string y se creo la nueva configuraci贸n con String solver
    public void generateReferencePoint(Organs[] o, int opt, String jobThreadID,String scoreConfigFile) throws IOException {

        switch (opt) {
            case 1:
                //this.solveLogFunction(o, jobThreadID,scoreConfigFile);
                break;
            case 2:
                this.solveLexiRectum(o, jobThreadID);
                break;
            case 4:
                this.solveVariableWeights(o, jobThreadID);
                break;
            case 5:
                this.solveConvexLogFunction(o, jobThreadID);
                break;
            case 6:
                //this.solveWeightedSum(o, jobThreadID,scoreConfigFile);//no existe, no confundir con solveWeightedSum(o, solver, jobThreadID)
                break;
            case 7:
                //this.solveAdaptiveWeightedSum(o, jobThreadID,scoreConfigFile);
                break;
            case 8:
                //this.solveULogFunction(o, jobThreadID, scoreConfigFile);
                break;
        }

    }
*/
    public void updateWeights(Organs[] o){//katty
        for (int i = 0; i < o.length; i++) {
            this.weights[i]=o[i].weight;
        }
    }
    
    public void calculateDoseAll(Organs[] o, String pathFile)
            throws IOException {//katty
        String dir;
        for (int i = 0; i < o.length; i++) {
            dir = pathFile + "DVH_" + o[i].name + ".txt";
            calculateDoseVolume(o,pathFile,dir, i, o[i].isTarget);
        }
    }
   
    public void calculateDose(Organs[] o,String dir, int i_organ, boolean isTarget)
            throws IOException {//katty {suma Aij*xj}d
        ArrayList<Double> testList = new ArrayList<Double>(); //dosis
        ArrayList<Double> volume = new ArrayList<Double>();  //volume
        //String sp="\t";
        //String dir = "./test2-OptDVH_Rectum.txt";
        File f = new File(dir);
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line = fileIn.readLine();
        int aux = 1;
        while (line != null) {
            if (!line.contains("%") && !line.contains("sum")) {
                //System.out.println(line.trim());
                String[] auxReader = line.split("\\s");
                for (int i = 0; i < auxReader.length; i++) {
                    if (!auxReader[i].equals("")) {
                        if ((aux % 2) == 0) {
                            //System.out.println("i= "+i+" value= "+auxReader[i]);
                            testList.add(Double.parseDouble(auxReader[i]));
                        }
                        aux++;
                    }
                }
            }
            line = fileIn.readLine();
        }
        fileIn.close();

        int testListSize = testList.size();//nro de voxels
        /*Double sizeAux = (double) testListSize;
        Double porcentaje = (double) 0;
        for(int i=0;i<testList.size();i++){
            porcentaje=((sizeAux-i)/sizeAux)*100;
            volume.add(porcentaje);
        }*/
        if (isTarget) {
            double D95 = 0;
            for (int i = 0; i < testListSize; i++) {
                if (testList.get(i) > (o[i_organ].desiredDose * 0.95)) {
                    D95=D95+1;
                }
            }
            this.Dxx[i_organ] = new double[1];
            D95 = (D95 / testListSize) * 100;
            this.Dxx[i_organ][0] = D95;
            System.out.println(o[i_organ].name+" D95 = "+D95);
        } else {
            double D75 = 0, D100 = 0, D50 = 0, Dmax = 0, Dmean = 0;
            for (int i = 0; i < testListSize; i++) {
                Dmean = Dmean + testList.get(i);
                if (testList.get(i) > Dmax) {
                    Dmax = testList.get(i);
                }
                if (testList.get(i) > (o[i_organ].desiredDose * 0.75)) {
                    D75++;
                }
                if (testList.get(i) > (o[i_organ].desiredDose * 0.50)) {
                    D50++;
                }
                if (testList.get(i) > (o[i_organ].desiredDose)) {
                    D100++;
                }
            }
            this.Dxx[i_organ] = new double[5];
            Dmean = Dmean / testListSize;
            D75 = (D75 / testListSize) * 100;
            D50 = (D50 / testListSize) * 100;
            D100 = (D100 / testListSize) * 100;
            this.Dxx[i_organ][0] = Dmean;
            this.Dxx[i_organ][1] = Dmax;
            this.Dxx[i_organ][2] = D50;
            this.Dxx[i_organ][3] = D75;
            this.Dxx[i_organ][4] = D100;
            System.out.println(o[i_organ].name+" Dmean = "+Dmean+" Dmax = "+Dmax+" D50 = "+D50+" D75 = "+D75+" D100 = "+D100);
        }

        //System.out.println(testList.get(i));
        //System.out.println("max= "+Dmax+" mean= "+Dmean+" 95= "+D95+" 75= "+D75+" 50= "+D50+" 100= "+D100);
    }
    
    public void calculateDoseVolume(Organs[] o,String jobThreadID, String dir, int i_organ, boolean isTarget)
            throws IOException {//katty {suma Aij*xj}d
        ArrayList<Double> testList = new ArrayList<>(); //dosis
        ArrayList<Double> volume = new ArrayList<>();  //volume
        //String sp="\t";
        //String dir = "./test2-OptDVH_Rectum.txt";
        //testList.clear() ;   
        //volume.clear();
        //testList=readFileDVH(dir) ;
        
        
        File f = new File(dir);
        
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line = fileIn.readLine();
        int aux = 1;
        while (line != null) {
            if (!line.contains("%") && !line.contains("sum")) {
                //System.out.println(line.trim());
                String[] auxReader = line.split("\\s");
                for (int i = 0; i < auxReader.length; i++) {
                    if (!auxReader[i].equals("")) {
                        if ((aux % 2) == 0) {
                            //System.out.println("i= "+i+" value= "+auxReader[i]);
                            testList.add(Double.parseDouble(auxReader[i]));
                           // System.out.println(Double.parseDouble(auxReader[i]));
                        }
                        aux++;
                    }
                }
            }
            line = fileIn.readLine();
        }
        fileIn.close();
        
        Collections.sort(testList);
        int testListSize = testList.size();//nro de voxels
        /**C贸digo para imprimir el DVH con Dosis y volumen*/
        Double porcentaje = (double) 0;
        Double sizeAux = (double) testListSize;
        double Dmean = 0;
        /*
        
        
        String newdirFile = jobThreadID + "DVH2_" + o[i_organ].name + ".txt";
        BufferedWriter bwFile = null;
        File solFile = new File(newdirFile);
        if (solFile.exists()) {
            if(!solFile.delete()){
                    System.out.println("Delete operation failed.(DVH2)");
            }
            //bwFile = new BufferedWriter(new FileWriter(solFile, true));
            bwFile = new BufferedWriter(new FileWriter(solFile));
        } else {
            bwFile = new BufferedWriter(new FileWriter(solFile));
        }
       // System.out.println("");
        //System.out.println("asdfffff");
        //System.out.println(testList.size()+"\t"+"\t"+o[i_organ].totalVoxels);
        for(int i=0;i<testList.size();i++){
            porcentaje=((sizeAux-i)/sizeAux)*100;
            volume.add(porcentaje);
           // System.out.println(testList.get(i) + "\t" + porcentaje );
            writeLine(testList.get(i) + "\t" + porcentaje + "\n", bwFile);
            
            //System.out.println(testList.get(i));
        }
        bwFile.close();
        */
        for(int i=0;i<testList.size();i++){
            porcentaje=((sizeAux-i)/sizeAux)*100;
            volume.add(porcentaje);
            Dmean = Dmean + testList.get(i);
           // System.out.println(testList.get(i) + "\t" + porcentaje );
           //  writeLine(testList.get(i) + "\t" + porcentaje + "\n", bwFile);
            
            //System.out.println(testList.get(i));
        }
        /**Fin elaboracion de archivo dvh**/
        if (isTarget) {
            double D95 = 0, V95 = 0,D100=0;
            int contador = 0;
            for (int i = 0; i < testListSize; i++) {
                double auxiliar=0;
                auxiliar=o[i_organ].desiredDose * 0.95;
                //System.out.println("auxiliar "+auxiliar);
               // if (testList.get(i) > auxiliar) {
               //     D95=D95+1;
               // }
                //if (testList.get(i) > o[i_organ].desiredDose) {
                //    D100=D100+1;
                //}
                if (testList.get(i) > auxiliar) {//https://www.ncbi.nlm.nih.gov/pubmed/16369743
                //    D95=D95+1;
                    V95=V95+volume.get(i);
                    contador++;
                    //System.out.println("contador: "+contador+" volumen "+volume.get(i)+" suma "+V95);
                }
            }
            this.Dxx[i_organ] = new double[7];
            //System.out.println("contador: "+contador+" suma volumen "+V95);
            if(contador>0){
               V95 =  V95 / contador ; 
            }else{
               V95 = 0; 
            }
            
            //D95 = ( D95 / testListSize ) * 100;
            //D100 = ( D100 / testListSize ) * 100;
            this.Dxx[i_organ][0] = getDxx_ge(testList,volume,95);//D95;
            this.Dxx[i_organ][1] = getDxx_ge(testList,volume,100);//D100;
            this.Dxx[i_organ][2] = V95;
            //String dir = pathFile + "DVH_" + o[o_index].name + ".txt";
            //ArrayList<Double> dose = readFileDVH(dir), volume = getVolume(dose);
            double D5, BED,b;
            D5 = getDxx_ge(testList,volume,5);
            D95 = getDxx_ge(testList,volume,95);
            Dmean = Dmean / testListSize;
            b=Dmean+Dmean/10;
            BED=b/(1+75/10);
            this.Dxx[i_organ][3] = Math.abs((o[i_organ].desiredDose-D95)+(D5-o[i_organ].desiredDose));
            this.Dxx[i_organ][4] = D5/D95;
            this.Dxx[i_organ][5] = D5;
            this.Dxx[i_organ][6] = BED;
            
            //for (int k=0; k<this.Dxx[i_organ].length;k++){
            //    this.Dxx[i_organ][k]=(double) (Math.floor(this.Dxx[i_organ][k] * 1000) / 1000);
            //}
//System.out.println(o[i_organ].name+" D95 = "+D95+" D100 = "+D100+" V95 = "+V95);
            System.out.println(o[i_organ].name+" D95 = "+this.Dxx[i_organ][0]+" D5 = "+this.Dxx[i_organ][5]+
                    " D100 = "+this.Dxx[i_organ][1]+" V95 = "+this.Dxx[i_organ][2]+" gEUDbased = "+this.Dxx[i_organ][3]+" H_index = "+this.Dxx[i_organ][4]+" BED = "+this.Dxx[i_organ][6]);
        } else {
            double D75 = 0, D100 = 0, D50 = 0, Dmax = 0;
            Dmean = 0;
            for (int i = 0; i < testListSize; i++) {
                Dmean = Dmean + testList.get(i);
                if (testList.get(i) > Dmax) {
                    Dmax = testList.get(i);
                }
                //if (testList.get(i) > (o[i_organ].desiredDose * 0.75)) {
                //    D75++;
                //}
                //if (testList.get(i) > (o[i_organ].desiredDose * 0.50)) {
                //    D50++;
                //}
                //if (testList.get(i) > (o[i_organ].desiredDose)) {
                //    D100++;
                //}
            }
            this.Dxx[i_organ] = new double[7];
            Dmean = Dmean / testListSize;
            //D75 = (D75 / testListSize) * 100;
            //D50 = (D50 / testListSize) * 100;
            //D100 = (D100 / testListSize) * 100;
            this.Dxx[i_organ][0] = Dmean;
            this.Dxx[i_organ][1] = Dmax;
            this.Dxx[i_organ][2] = getDxx_ge(testList,volume,50);//D50;
            this.Dxx[i_organ][3] = getDxx_ge(testList,volume,75);//D75;
            this.Dxx[i_organ][4] = getDxx_ge(testList,volume,100);//D100;
            this.Dxx[i_organ][5] = Dmax/o[i_organ].doseUB;
            this.Dxx[i_organ][6] = Dmean/o[i_organ].desiredDose;
            
            //for (int k=0; k<this.Dxx[i_organ].length;k++){
            //    this.Dxx[i_organ][k]=(double) (Math.floor(this.Dxx[i_organ][k] * 1000) / 1000);
            //}
            
            //System.out.println(o[i_organ].name+" Dmean = "+Dmean+" Dmax = "+Dmax+" D50 = "+D50+" D75 = "+D75+" D100 = "+D100);
            System.out.println(o[i_organ].name+" Dmean = "+this.Dxx[i_organ][0]
                +" Dmax = "+this.Dxx[i_organ][1]+" D50 = "+this.Dxx[i_organ][2]+" D75 = "+this.Dxx[i_organ][3]+" D100 = "+this.Dxx[i_organ][4]+" Dmax/Dc = "+this.Dxx[i_organ][5]+" Dmean/Dc = "+this.Dxx[i_organ][6]);
     
        }
        testList.clear();
        volume.clear();
        //System.out.println(testList.get(i));
        //System.out.println("max= "+Dmax+" mean= "+Dmean+" 95= "+D95+" 75= "+D75+" 50= "+D50+" 100= "+D100);
    }
   
    public double calculateLTCP(double alpha, Organs[] o,String pathFile,int o_index)throws IOException{//katy revisar
        double LTCP=0;
        //m=o[i].totalVoxels
        //PTV  valores  input file agregar el alpha y guardarlo
        Double Do;
        Integer m=0;
        String dir;
        ArrayList<Double> dvhList = new ArrayList<Double>(); //dosis
        
        
       // for(int i=0;i<o.length;i++){
            if(o[o_index].isTarget){
                m=o[o_index].totalVoxels;
                dir = pathFile + "DVH_" + o[o_index].name + ".txt";
                File f = new File(dir);
                BufferedReader fileIn = new BufferedReader(new FileReader(f));
                String line = "";
                line = fileIn.readLine();
                int aux = 1;
                while (line != null) {
                    if (!line.contains("%") && !line.contains("sum")) {
                        //System.out.println(line.trim());
                        String[] auxReader = line.split("\\s");
                        for (int k = 0; k < auxReader.length; k++) {
                            if (!auxReader[k].equals("")) {
                                if ((aux % 2) == 0) {
                                    //System.out.println("i= "+i+" value= "+auxReader[i]);
                                    dvhList.add(Double.parseDouble(auxReader[k]));
                                }
                                aux++;
                            }
                        }
                    }
                    line = fileIn.readLine();
                }
                fileIn.close();
                //Do = calculateObtainedDose(o,o_index);
                Do=o[o_index].desiredDose;
                for(int j=0;j<m;j++){
                    LTCP=LTCP+Math.exp(-alpha*(dvhList.get(j)-Do));
                }
            }    
       // }
        LTCP=LTCP/m;
        return (LTCP);
    }
    
    public double calculateTCP(double TCD50,double D,double gama50){
        double TCP=0;
        TCP=Math.pow(1+Math.pow((TCD50/D),4*gama50), -1);
        return (TCP);
    }
    
    public  double calculateNTCP(double TD50,double D,double gama50){
        double NTCP;
        NTCP=Math.pow(1+Math.pow((TD50/D),4*gama50), -1);
        return (NTCP);
    }
  
    public void calculateAllScores(Organs[] o,String pathFile
            ,String [][] scoreNames,double [][][] scoreConfig,int n_targetScores
            ,int n_oarScores ) throws IOException{
        Double Do;
        for(int o_index = 0 ; o_index < o.length; o_index++) {
            Do = calculateObtainedDose(o,o_index);
            if(o[o_index].isTarget){
                this.allScores[o_index]=new double [n_targetScores];
                for(int j = 0; j < n_targetScores ; j++){
                    if(scoreNames[o_index][j].equals("LTCP")){
                        this.allScores[o_index][j] = calculateLTCP(scoreConfig[o_index][j][0],
                                o,pathFile,o_index);
                        //this.allScores[o_index][j]=(double) (Math.floor(this.allScores[o_index][j]* 1000) / 1000);
                    }
                    if(scoreNames[o_index][j].equals("TCP")){
                        this.allScores[o_index][j] = calculateTCP(scoreConfig[o_index][j][1],
                            this.gEUD[o_index],scoreConfig[o_index][j][0]);
                        //this.allScores[o_index][j]=(double) (Math.floor(this.allScores[o_index][j]* 1000) / 1000);
                    }
                    if(scoreNames[o_index][j].equals("default")){
                        this.allScores[o_index][j] = (o[o_index].desiredDose / Do);
                      //  this.allScores[o_index][j]=(double) (Math.floor(this.allScores[o_index][j]* 1000) / 1000);
                    }
                    System.out.print(o[o_index].name+" "+scoreNames[o_index][j]+": "+this.allScores[o_index][j]+" ");
                    //[2] default en base al paper globalScore de Rocha y Dias
                    /*
                    for(int i=0;i<targetScores;i++){
                        System.out.print(scoreNames[o_index][i]+": "+this.allScores[o_index][i]);
                    }*/
                }
            }else{
                this.allScores[o_index] = new double [n_oarScores];
                
                for(int j = 0; j < n_oarScores ; j++){
                    if(scoreNames[o_index][j].equals("NTCP")){
                        this.allScores[o_index][j] = calculateNTCP(scoreConfig[o_index][j][1],
                        this.gEUD[o_index],scoreConfig[o_index][j][0]);
                       // this.allScores[o_index][j]=(double) (Math.floor(this.allScores[o_index][j]* 1000) / 1000);
                    }
                    if(scoreNames[o_index][j].equals("default")){
                        this.allScores[o_index][j] = Do / o[o_index].desiredDose; //default
                       // this.allScores[o_index][j]=(double) (Math.floor(this.allScores[o_index][j]* 1000) / 1000);
                    }
                    System.out.print(o[o_index].name+" "+scoreNames[o_index][j]+": "+this.allScores[o_index][j]+" ");
                    /*for(int i=0;i<oarScores;i++){
                        System.out.print(scoreNames[o_index][i]+": "+this.allScores[o_index][i]);
                    }*/
                }
            }
            System.out.println(" ");
        }
    }
    public static double getDxx_ge(ArrayList<Double> dose,ArrayList<Double> volume,double percentage){
	double Dxx = (double)-1;
	int i=0;
	for(i=0;i<volume.size();i++){
	//    System.out.println("dose: "+dose.get(i)+"\t volume: "+volume.get(i)+">="+percentage);
		if(volume.get(i)<percentage){
			return dose.get(i-1);
		}
	}
	if(percentage==0){
		return dose.get(i-1);
	}
	return Dxx;
}
public static double getDxx_le(ArrayList<Double> dose,ArrayList<Double> volume,double percentage){
	double Dxx = (double)-1;

	int i=0;
	for(i=0;i<volume.size();i++){
	//    System.out.println("dose: "+dose.get(i)+"\t volume: "+volume.get(i)+"<="+percentage);
		if(volume.get(i)<=percentage){
			return dose.get(i);
		}
	}
	if(percentage==0){
		return dose.get(i-1);
	}
        return Dxx;
}
    public static ArrayList<Double> readFileDVH(String dir) 
            throws FileNotFoundException, IOException{
        ArrayList<Double> testList = new ArrayList<Double>(); //dosis
        //String sp="\t";
        //String dir = "./test2-OptDVH_Rectum.txt";
        File f = new File(dir);
        try (BufferedReader fileIn = new BufferedReader(new FileReader(f))) {
            String line = "";
            line = fileIn.readLine();
            int aux = 1;
            while (line != null) {
                if (!line.contains("%") && !line.contains("sum")) {
                    String[] auxReader = line.split("\\s");
                    for (int i = 0; i < auxReader.length; i++) {
                        if (!auxReader[i].equals("")) {
                            if ((aux % 2) == 0) {
                                //System.out.println("i= "+i+" value= "+auxReader[i]);
                                testList.add(Double.parseDouble(auxReader[i]));
                            }
                            aux++;
                        }
                    }
                }
                line = fileIn.readLine();
            }
            fileIn.close();
        }
        Collections.sort(testList);
        return testList;
    }
    public ArrayList<Double> getVolume(ArrayList<Double> dose){
        int doseSize = dose.size();//nro de voxels
        ArrayList<Double> volume = new ArrayList<>();  //volume
        /**C贸digo para imprimir el DVH con Dosis y volumen*/
        Double sizeAux = (double) doseSize;
        Double porcentaje = (double) 0;
        /*String newdirFile = jobThreadID + "DVH2_" + o[i_organ].name + ".txt";
        BufferedWriter bwFile = null;
        File solFile = new File(newdirFile);
        if (solFile.exists()) {
            if(!solFile.delete()){
                    System.out.println("Delete operation failed.(DVH2)");
            }
            //bwFile = new BufferedWriter(new FileWriter(solFile, true));
            bwFile = new BufferedWriter(new FileWriter(solFile));
        } else {
            bwFile = new BufferedWriter(new FileWriter(solFile));
        }*/
       // System.out.println("");
        //System.out.println("asdfffff");
        //System.out.println(dose.size()+"\t"+"\t"+o[i_organ].totalVoxels);
        for(int i=0;i<dose.size();i++){
            porcentaje=((sizeAux-i)/sizeAux)*100;
            volume.add(porcentaje);
           // System.out.println(dose.get(i) + "\t" + porcentaje );
           // writeLine(dose.get(i) + "\t" + porcentaje + "\n", bwFile);
            
            //System.out.println(dose.get(i));
        }
        //bwFile.close();
        /**Fin elaboracion de archivo dvh**/
        return volume;
    }
    
    public void calculateSingleScore(Organs[] o,String pathFile,String scoreType
    ,ArrayList<Double> parametersList,int o_index) throws IOException {//katty
        Double Do;
        //Double aux;
        Do = calculateObtainedDose(o,o_index);
        if(o[o_index].isTarget){
            //System.out.println("LTCP: "+calculateLTCP(parametersList.get(0),o,pathFile,o_index));
            //System.out.println("TCP: "+calculateTCP(parametersList.get(1),this.gEUD[o_index],parametersList.get(0)));
            //aux=(o[o_index].desiredDose / Do);
            //System.out.println("Target Score: "+aux);
            switch(scoreType){
                case "LTCP":
                    this.scores[o_index] = calculateLTCP(parametersList.get(0),o,pathFile,o_index);
                    System.out.println("LTCP: "+this.scores[o_index]);  
                    break;
                case "TCP":
                    this.scores[o_index] = calculateTCP(parametersList.get(1),this.gEUD[o_index],parametersList.get(0));
                    System.out.println("TCP: "+this.scores[o_index]);
                    break;
                case "default":
                    this.scores[o_index] = o[o_index].desiredDose / Do;
                    System.out.println("Target Score: "+this.scores[o_index]);
                    break;
                case "D95":
                    this.scores[o_index]=Dxx[o_index][0];
                    System.out.println("D95: "+this.scores[o_index]);
                    break;
                case "D100":
                    this.scores[o_index]=Dxx[o_index][1];
                    System.out.println("D100: "+this.scores[o_index]);
                    break;
                case "V95":
                    this.scores[o_index]=Dxx[o_index][2];
                    System.out.println("V95: "+this.scores[o_index]);
                    break;
                case "gEUDbased"://Math.abs((o[i_organ].desiredDose-D95)+(D5-o[i_organ].desiredDose));
                    this.scores[o_index]=this.Dxx[o_index][3];
                    System.out.println("gEUDbased: "+this.scores[o_index]);
                    break;
                case "EUDbased"://Math.abs((o[i_organ].desiredDose-D95)+(D5-o[i_organ].desiredDose));
                    this.scores[o_index]=this.Dxx[o_index][3];
                    System.out.println("EUDbased: "+this.scores[o_index]);
                    break;
                case "H_index"://D5/D95
                    this.scores[o_index]=this.Dxx[o_index][4];
                    System.out.println("H_index: "+this.scores[o_index]);
                    break;
                case "BED"://D5/D95
                    this.scores[o_index]=this.Dxx[o_index][5];
                    System.out.println("BED: "+this.scores[o_index]);
                    break;
                default:
                    this.scores[o_index] = o[o_index].desiredDose / Do;
                    System.out.println("Target Score: "+this.scores[o_index]);
                    break;   
            }/*
            if(scoreType.equals("LTCP")){
                this.scores[o_index] = calculateLTCP(parametersList.get(0),o,pathFile,o_index);
                System.out.println("LTCP: "+this.scores[o_index]);
            }else{
                if(scoreType.equals("TCP")){
                    this.scores[o_index] = calculateTCP(parametersList.get(1),this.gEUD[o_index],parametersList.get(0));
                    System.out.println("TCP: "+this.scores[o_index]);
                }else{
                    if(scoreType.equals("default")){
                        this.scores[o_index] = o[o_index].desiredDose / Do;
                        System.out.println("Target Score: "+this.scores[o_index]);
                    }
                }
            }*/
        }else{
            //System.out.println("NTCP: "+calculateNTCP(parametersList.get(1),this.gEUD[o_index],parametersList.get(0)));
            //aux=Do / o[o_index].desiredDose;
            //System.out.println("OAR Score: "+aux);
            switch(scoreType){
                case "NTCP":
                    this.scores[o_index] = calculateNTCP(parametersList.get(1),this.gEUD[o_index],parametersList.get(0));
                    System.out.println("NTCP: "+this.scores[o_index]);
                    break;
                case "default":
                    this.scores[o_index] = Do / o[o_index].desiredDose;
                    System.out.println("OAR Score: "+this.scores[o_index]);
                    break;
                case "Dmean":
                    this.scores[o_index]=Dxx[o_index][0];
                    System.out.println("Dmean: "+this.scores[o_index]);
                    break;
                case "Dmax":
                    this.scores[o_index]=Dxx[o_index][1];
                    System.out.println("Dmax: "+this.scores[o_index]);
                    break;
                case "D50":
                    this.scores[o_index]=Dxx[o_index][2];
                    System.out.println("D50: "+this.scores[o_index]);
                    break;
                case "D75":
                    this.scores[o_index]=Dxx[o_index][3];
                    System.out.println("D75: "+this.scores[o_index]);
                    break;
                case "D100":
                    this.scores[o_index]=Dxx[o_index][4];
                    System.out.println("D100: "+this.scores[o_index]);
                    break;
                case "DmaxDc":
                    this.scores[o_index]=Dxx[o_index][5];
                    System.out.println("Dmax/Dc: "+this.scores[o_index]);
                    break;
                case "DmeanDc":
                    this.scores[o_index]=Dxx[o_index][6];
                    System.out.println("Dmean/Dc: "+this.scores[o_index]);
                    break;
                default:
                    this.scores[o_index] = Do / o[o_index].desiredDose;
                    System.out.println("OAR Score: "+this.scores[o_index]);
                    break;
            }/*
            if(scoreType.equals("NTCP")){
                this.scores[o_index] = calculateNTCP(parametersList.get(1),this.gEUD[o_index],parametersList.get(0));
                System.out.println("NTCP: "+this.scores[o_index]);
            }else{
                if(scoreType.equals("default")){
                    this.scores[o_index] = Do / o[o_index].desiredDose;
                    System.out.println("OAR Score: "+this.scores[o_index]);
                }
            }*/
        }
    }
	
    public double calculateObtainedDose(Organs[] o,int i_organ) throws IOException {//katty
        double Dc = 0;
        if(o[i_organ].isTarget){
            Dc = this.Dxx[i_organ][0];//D95
        }else{
            //use max or mean tolerance doses 
            //depending on the type of organ(serial or parallel)
            //Dc = this.Dxx[i_organ][0];//mean
            //Dc = this.Dxx[i_organ][1];//max
            if(o[i_organ].name.equals("RECTUM")||o[i_organ].name.equals("RECTO")){//serial
                Dc = this.Dxx[i_organ][1];//max
            }else{// If it is Bladder is parallel
                Dc = this.Dxx[i_organ][0];//mean    
            }
        }
        return Dc;
    }
    
    public void calculateGlobalScore(Organs[] o,ArrayList<Double> weightsList) throws IOException {//katty
        double aux_globalScore = 0;
        for (int i = 0; i < o.length; i++) {
            aux_globalScore = aux_globalScore + (this.scores[i] * weightsList.get(i));
        }
        this.globalScore = aux_globalScore;
        System.out.println("Global Score = " + this.globalScore);
    }
    
    public void solveLogFunction(Organs[] o, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];
        updateWeights(o);//katty
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver logFunction = new AMPL_Solver(o, this, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            logFunction.generateParametersFile_logFunction(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            logFunction.generateParametersFile_logFunction();
        }
        logFunction.runLogFunction_Solver(o);
        logFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(logFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getLogFunctionValue(o);
    }
    
    public void solveLogFunction(Organs[] o, String solver, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];
        updateWeights(o);//katty
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver logFunction = new AMPL_Solver(o, this, solver, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            logFunction.generateParametersFile_logFunction(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            logFunction.generateParametersFile_logFunction();
        }
        logFunction.runLogFunction_Solver(o);
        logFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(logFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getLogFunctionValue(o);
    }


    //public void solveLogFunction(Organs_bkp[] o, DDM_bkp M) throws IOException{
    public void solveLogFunction(Organs[] o, String solver, String jobThreadID,String scoreConfigFile) 
            throws IOException {
        double[] EUD0 = new double[o.length];
        //updateWeights(o);//katty
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver logFunction = new AMPL_Solver(o, this,solver, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            logFunction.generateParametersFile_logFunction(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            logFunction.generateParametersFile_logFunction();
        }
        logFunction.runLogFunction_Solver(o);
        logFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(logFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getLogFunctionValue(o);
        //agregar llamadas a funciones katty
        this.calculateDoseAll(o,jobThreadID);//katy
        this.readScoresConfig(scoreConfigFile,jobThreadID,o);
        //this.calculateSingleScore(o,jobThreadID);//katy
        //this.calculateGlobalScore(o);//katy
        
    }
    
    public void solveULogFunction(Organs[] o, String solver, String jobThreadID) 
            throws IOException {
        double[] EUD0 = new double[o.length];
        //updateWeights(o);//katty
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver UlogFunction = new AMPL_Solver(o, this,solver, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            UlogFunction.generateParametersFile_UlogFunction(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            UlogFunction.generateParametersFile_UlogFunction();
        }
        UlogFunction.runULogFunction_Solver(o);
        UlogFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(UlogFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getULogFunctionValue(o);
        //agregar llamadas a funciones katty
        //this.calculateDoseAll(o,jobThreadID);//katy
        //this.readScoresConfig(scoreConfigFile,jobThreadID,o);
        //this.calculateSingleScore(o,jobThreadID);//katy
        //this.calculateGlobalScore(o);//katy
        
    }
    
    public void solveULogFunction(Organs[] o, String solver, String jobThreadID,String scoreConfigFile) 
            throws IOException {
        double[] EUD0 = new double[o.length];
        //updateWeights(o);//katty
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver UlogFunction = new AMPL_Solver(o, this,solver, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            UlogFunction.generateParametersFile_UlogFunction(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            UlogFunction.generateParametersFile_UlogFunction();
        }
        UlogFunction.runULogFunction_Solver(o);
        UlogFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(UlogFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getULogFunctionValue(o);
        //agregar llamadas a funciones katty
        this.calculateDoseAll(o,jobThreadID);//katy
        this.readScoresConfig(scoreConfigFile,jobThreadID,o);
        //this.calculateSingleScore(o,jobThreadID);//katy
        //this.calculateGlobalScore(o);//katy
        
    }
    public void solveCuadraticFunction(Organs[] o, String solver, double maxIntensity , String jobThreadID, int weighted) throws IOException { //Maicholl
        double[] EUD0 = new double[o.length];
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }

        Random r = new Random();
        this.weights[1] = 0;
        while (this.weights[1] == 0 || this.weights[1] == 1) {
            this.weights[1] = Math.floor(100.0 * r.nextDouble()) / 100.0;
        }
        this.weights[2] = 1 - this.weights[1];

        AMPL_Solver cuadraticFunction = new AMPL_Solver(o, this, solver, maxIntensity, 0, jobThreadID);
        cuadraticFunction.wPar[1] = this.weights[1];
        cuadraticFunction.wPar[2] = this.weights[2];
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            cuadraticFunction.generateParametersFile_cuadraticSum(this.intensity);
            System.out.println("Generating parameters Files for JobID "+jobThreadID+" (using previous intensities)....Done!");
        } else {
            cuadraticFunction.generateParametersFile_cuadraticSum();
            System.out.println("Generating parameters Files for JobID "+jobThreadID+" (without using previous intensities)...Done!");
        }
        //Delete solution files (if any) from previous iteration
        cuadraticFunction.deleteSolutionFiles(o);
        
        
        cuadraticFunction.runCuadraticSum_Solver(o,weighted);
        
        cuadraticFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(cuadraticFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getCuadraticSumValue(jobThreadID,o);
        
    }
    
    public void solveCuadraticFunction(Organs[] o, String solver, double maxIntensity , String jobThreadID,String scoreConfigFile) throws IOException { //Maicholl
        double[] EUD0 = new double[o.length];
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }

        Random r = new Random();
        this.weights[1] = 0;
        while (this.weights[1] == 0 || this.weights[1] == 1) {
            this.weights[1] = Math.floor(100.0 * r.nextDouble()) / 100.0;
        }
        this.weights[2] = 1 - this.weights[1];
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver cuadraticFunction = new AMPL_Solver(o, this, solver, maxIntensity, 0, jobThreadID);
        cuadraticFunction.wPar[1] = this.weights[1];
        cuadraticFunction.wPar[2] = this.weights[2];
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            cuadraticFunction.generateParametersFile_cuadraticSum(this.intensity);
            System.out.println("Generating parameters Files for JobID "+jobThreadID+" (using previous intensities)....Done!");
        } else {
            cuadraticFunction.generateParametersFile_cuadraticSum();
            System.out.println("Generating parameters Files for JobID "+jobThreadID+" (without using previous intensities)...Done!");
        }
        //Delete solution files (if any) from previous iteration
        cuadraticFunction.deleteSolutionFiles(o);
        
        cuadraticFunction.runCuadraticSum_Solver(o);
        cuadraticFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(cuadraticFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getCuadraticSumValue(jobThreadID,o);
        //agregar llamadas a funciones katty
        this.calculateDoseAll(o,jobThreadID);//katy
        this.readScoresConfig(scoreConfigFile,jobThreadID,o);
    }
    public void readScoresConfig(String dir,String jobThreadID,Organs[] o) throws IOException{
        Integer numOrgans = null;
        String sp="\t";
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
        Integer total=0;    
        String scoreType;
        Integer o_index;
        ArrayList<Double> parametersList=new ArrayList<>();
        ArrayList<Double> weightsList=new ArrayList<>();
        while(line != null){
            if (!line.contains("%")){
                for (int y=0;y<numOrgans;y++){
                    String[] auxReader = line.split(sp);
                    o_index=Integer.parseInt(auxReader[0]);
                    scoreType=auxReader[1];
                    //weight falta
                    /*System.out.println("index "+Integer.parseInt(auxReader[0])+
                            " Score-type "+ auxReader[1]+" weight "+
                            Double.parseDouble(auxReader[2]));*/
                    weightsList.add(Double.parseDouble(auxReader[2]));
                    total=0;
                    if( auxReader[1].equals("LTCP")){
                        total=1+3;
                    }else{ 
                        if(auxReader[1].equals("TCP") || auxReader[1].equals("NTCP")){
                            total=2+3;
                        }
                    }
                    for(int i=3;i<total;i++){
                        //System.out.println("Parameter "+Double.parseDouble(auxReader[i]));
                        parametersList.add(Double.parseDouble(auxReader[i]));
                    }
                    
                    this.calculateSingleScore(o,jobThreadID,scoreType,parametersList,o_index);
                    parametersList.clear();
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        //#scores para target y #scores oar
        int n_targetScores = 0;
        int n_oarScores = 0;
        int n_Parameters = 0;
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                n_targetScores = Integer.parseInt(auxReader[0]);
                n_oarScores = Integer.parseInt(auxReader[1]);
                n_Parameters = Integer.parseInt(auxReader[2]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //readAllScoresConfig
        String scoreNames[][] = new String [o.length][];
        double scoreConfig[][][] = new double [o.length][][];
        int aux;
        aux = -1;
        while(line != null){
            if (!line.contains("%")){
                //%i	j	name	#parameters parameters
                String[] auxReader = line.split(sp);
                
                int i = Integer.parseInt(auxReader[0]);//i
                int j = Integer.parseInt(auxReader[1]);//j
                
                //Double.parseDouble(auxReader[2]);//scoreName
                int nParameters = Integer.parseInt(auxReader[3]);//#parameters
                if(aux!=i){
                    if(o[i].isTarget){
                        aux = i;
                        scoreNames[i] = new String [n_targetScores];
                        scoreConfig[i] = new double [n_targetScores][n_Parameters];
                    }else{
                        aux = i;
                        scoreNames[i] = new String [n_oarScores];
                        scoreConfig[i] = new double [n_oarScores][n_Parameters];
                    }
                }
               
                scoreNames[i][j] = auxReader[2];
                nParameters=4+nParameters;
                int index=0;
                for (int k=4;k<nParameters;k++){
                    
                    scoreConfig[i][j][index]=Double.parseDouble(auxReader[k]);
                    index++;
                }
                //line=fileIn.readLine();
                //break;
            }
            line=fileIn.readLine();
        }
        fileIn.close();
        this.calculateAllScores(o, jobThreadID, scoreNames,scoreConfig,n_targetScores,n_oarScores);
        this.calculateGlobalScore(o,weightsList);
    }
    
   
    public void solveUnconstrainedWeightedSum_gEUD(Organs[] o, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        /**
         * ****************
         */
        /*Determine Weights*/
        /**
         * ****************
         */

        Random r = new Random();//quitar hacer update
        this.weights[1] = 0;
        while (this.weights[1] == 0 || this.weights[1] == 1) {
            this.weights[1] = Math.floor(100.0 * r.nextDouble()) / 100.0;
        }
        this.weights[2] = 1 - this.weights[1];
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver variableWeights = new AMPL_Solver(o, this, 0, jobThreadID);
        variableWeights.wPar[1] = this.weights[1];
        variableWeights.wPar[2] = this.weights[2];
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            variableWeights.generateParametersFile_unconstrainedWeightedSum_gEUD(this.intensity);//generateParametersFile_variableWeights(this.intensity);//Descomentar
        } else {
            System.out.println("Generating parameters File (1)");
            variableWeights.generateParametersFile_unconstrainedWeightedSum_gEUD();//generateParametersFile_variableWeights();//Descomentar
        }
        variableWeights.runVariableWeights_Solver(o);
        variableWeights.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(variableWeights.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getVariableWeightsValue(o);
    }
    
    
    public void solveUnconstrainedWeightedSum_gEUD(Organs[] o,String solver, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
       // /**
       //  * ****************
        // */
        ///*Determine Weights*/
        ///**
         //* ****************
         //*/

        /*Random r = new Random();
        this.weights[1] = 0;
        while (this.weights[1] == 0 || this.weights[1] == 1) {
            this.weights[1] = Math.floor(100.0 * r.nextDouble()) / 100.0;
        }
        this.weights[2] = 1 - this.weights[1];
        */
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver variableWeights = new AMPL_Solver(o, this, 0, jobThreadID);
        variableWeights.wPar[1] = this.weights[1];
        variableWeights.wPar[2] = this.weights[2];
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            variableWeights.generateParametersFile_unconstrainedWeightedSum_gEUD(this.intensity);//generateParametersFile_variableWeights(this.intensity);//descomentar
        } else {
            System.out.println("Generating parameters File (1)");
            variableWeights.generateParametersFile_unconstrainedWeightedSum_gEUD();//generateParametersFile_variableWeights();//descomentar
        }
        variableWeights.runUnconstrainedWeightedSum_gEUD_Solver(o);
        variableWeights.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(variableWeights.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getUnconstrainedWeightedSum_gEUD_Value(o);
    }
    
    public void solveUnconstrainedWeightedSum_gEUD(Organs[] o,String solver, String jobThreadID,String scoreConfigFile) throws IOException {
        double[] EUD0 = new double[o.length];

        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver UnconstrainedWeightedSum_gEUD = new AMPL_Solver(o, this,solver, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            UnconstrainedWeightedSum_gEUD.generateParametersFile_unconstrainedWeightedSum_gEUD(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            UnconstrainedWeightedSum_gEUD.generateParametersFile_unconstrainedWeightedSum_gEUD();
        }
        UnconstrainedWeightedSum_gEUD.runUnconstrainedWeightedSum_gEUD_Solver(o);//descomentar
        UnconstrainedWeightedSum_gEUD.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(UnconstrainedWeightedSum_gEUD.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getUnconstrainedWeightedSum_gEUD_Value(o);//descomentar
        this.calculateDoseAll(o,jobThreadID);//katy
        this.readScoresConfig(scoreConfigFile,jobThreadID,o);
    }

    public void solveAdaptiveWeightedSum(Organs[] o, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver adaptiveWeightedSum = new AMPL_Solver(o, this, 0, jobThreadID);
        adaptiveWeightedSum.wPar[1] = this.weights[1];
        adaptiveWeightedSum.wPar[2] = this.weights[2];
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            adaptiveWeightedSum.generateParametersFile_adaptiveWeightedSum(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            adaptiveWeightedSum.generateParametersFile_adaptiveWeightedSum();
        }
        adaptiveWeightedSum.runAdaptiveWeightedSum_Solver(o);
        adaptiveWeightedSum.getSolution();
        System.arraycopy(adaptiveWeightedSum.x, 0, this.intensity, 0, this.beamlets);
        //Random r = new Random();
        //this.intensity=new double[this.beamlets];
        //for(int i=0;i<this.intensity.length;i++){
        //    this.intensity[i] = r.nextDouble();
        //}

        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        //this.gEUD[0]=70;
        //this.gEUD[1] = 51.0 + r.nextDouble() * 8.0;
        //this.gEUD[2] = 31.0 + r.nextDouble() * 5.0;

        this.getAdaptiveWeightedSumValue(o);
    }
    
    public void solveAdaptiveWeightedSum(Organs[] o,String solver, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver adaptiveWeightedSum = new AMPL_Solver(o, this,solver, 0, jobThreadID);
        adaptiveWeightedSum.wPar[1] = this.weights[1];
        adaptiveWeightedSum.wPar[2] = this.weights[2];
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            adaptiveWeightedSum.generateParametersFile_adaptiveWeightedSum(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            adaptiveWeightedSum.generateParametersFile_adaptiveWeightedSum();
        }
        adaptiveWeightedSum.runAdaptiveWeightedSum_Solver(o);
        adaptiveWeightedSum.getSolution();
        System.arraycopy(adaptiveWeightedSum.x, 0, this.intensity, 0, this.beamlets);
        //Random r = new Random();
        //this.intensity=new double[this.beamlets];
        //for(int i=0;i<this.intensity.length;i++){
        //    this.intensity[i] = r.nextDouble();
        //}

        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        //this.gEUD[0]=70;
        //this.gEUD[1] = 51.0 + r.nextDouble() * 8.0;
        //this.gEUD[2] = 31.0 + r.nextDouble() * 5.0;

        this.getAdaptiveWeightedSumValue(o);
    }
    
    public void solveAdaptiveWeightedSum(Organs[] o, String solver,String jobThreadID,String scoreConfigFile) throws IOException {
        double[] EUD0 = new double[o.length];
        //updateWeights(o);//katty
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver adaptiveWeightedSum = new AMPL_Solver(o, this, 0, jobThreadID);
        adaptiveWeightedSum.wPar[1] = this.weights[1];
        adaptiveWeightedSum.wPar[2] = this.weights[2];
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            adaptiveWeightedSum.generateParametersFile_adaptiveWeightedSum(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            adaptiveWeightedSum.generateParametersFile_adaptiveWeightedSum();
        }
        adaptiveWeightedSum.runAdaptiveWeightedSum_Solver(o);
        adaptiveWeightedSum.getSolution();
        System.arraycopy(adaptiveWeightedSum.x, 0, this.intensity, 0, this.beamlets);
        //Random r = new Random();
        //this.intensity=new double[this.beamlets];
        //for(int i=0;i<this.intensity.length;i++){
        //    this.intensity[i] = r.nextDouble();
        //}

        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        //this.gEUD[0]=70;
        //this.gEUD[1] = 51.0 + r.nextDouble() * 8.0;
        //this.gEUD[2] = 31.0 + r.nextDouble() * 5.0;

        this.getAdaptiveWeightedSumValue(o);
        //agregar llamadas a funciones katty
        this.calculateDoseAll(o,jobThreadID);//katy
        this.readScoresConfig(scoreConfigFile,jobThreadID,o);
        //this.calculateSingleScore(o,jobThreadID);//katy
        //this.calculateGlobalScore(o);//katy
    }
    public void solveConvexLogFunction(Organs[] o, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];

        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver convexLogFunction = new AMPL_Solver(o, this, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            convexLogFunction.generateParametersFile_convexLogFunction(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            convexLogFunction.generateParametersFile_convexLogFunction();
        }
        convexLogFunction.runConvexLogFunction_Solver(o);
        convexLogFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(convexLogFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getConvexLogFunctionValue(o);
    }

    public void solveConvexLogFunction(Organs[] o, String solver,String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];

        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver convexLogFunction = new AMPL_Solver(o, this,solver, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            convexLogFunction.generateParametersFile_convexLogFunction(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            convexLogFunction.generateParametersFile_convexLogFunction();
        }
        convexLogFunction.runConvexLogFunction_Solver(o);
        convexLogFunction.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(convexLogFunction.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getConvexLogFunctionValue(o);
    }

    public void solveWeightedSum(Organs[] o, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];

        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver weightedSum = new AMPL_Solver(o, this, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            weightedSum.generateParametersFile_weightedSum(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            weightedSum.generateParametersFile_weightedSum();
        }
        weightedSum.runWeightedSum_Solver(o);
        weightedSum.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(weightedSum.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getWeightedSumValue(o);
    }

    public void solveWeightedSum(Organs[] o,String solver, String jobThreadID) throws IOException {
        double[] EUD0 = new double[o.length];

        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver weightedSum = new AMPL_Solver(o, this,solver, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            weightedSum.generateParametersFile_weightedSum(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            weightedSum.generateParametersFile_weightedSum();
        }
        weightedSum.runWeightedSum_Solver(o);
        weightedSum.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(weightedSum.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getWeightedSumValue(o);
    }

    public void solveTestFunction(Organs[] o, String jobThreadID) throws IOException {
        Random r = new Random();
        double aux;
        aux = 0;
        this.intensity = new double[this.beamlets];
        for (int i = 0; i < this.beamlets; i++) {
            this.intensity[i] = r.nextDouble();
        }
        this.gEUD[0] = o[0].desiredDose;
        aux = (1 + r.nextDouble());
        this.gEUD[1] = o[1].desiredDose * aux;
        aux = (1 + r.nextDouble());
        this.gEUD[2] = o[2].desiredDose * aux;
        aux = r.nextDouble();
        this.gEUD[3] = o[0].doseUB * aux;
        aux = 1 + r.nextDouble();
        this.singleObjectiveValue = aux * 5;
    }
    public void solveWeightedSum(Organs[] o,String solver, String jobThreadID,String scoreConfigFile) throws IOException {
        double[] EUD0 = new double[o.length];
       // updateWeights(o);//katty
        for (int i = 0; i < o.length; i++) {
            if (o[i].isTarget) {
                EUD0[i] = o[i].doseLB; //Target
            } else {
                EUD0[i] = o[i].doseUB; //OAR
            }
        }
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver weightedSum = new AMPL_Solver(o, this,solver, 0, jobThreadID);
        /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
        double auxSumIntensity = 0;
        for (int i = 0; i < this.intensity.length; i++) {
            auxSumIntensity = +this.intensity[i];
        }
        if (auxSumIntensity > 0) {
            System.out.println("Generating parameters File (0)");
            weightedSum.generateParametersFile_weightedSum(this.intensity);
        } else {
            System.out.println("Generating parameters File (1)");
            weightedSum.generateParametersFile_weightedSum();
        }
        weightedSum.runWeightedSum_Solver(o);
        weightedSum.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(weightedSum.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.getWeightedSumValue(o);
        //katy
        this.calculateDoseAll(o,jobThreadID);//katy
        this.readScoresConfig(scoreConfigFile,jobThreadID,o);
    }
    //public void getLogFunctionValue(Organs_bkp[] o) throws IOException{
    public void getLogFunctionValue(Organs[] o) throws IOException {
        double objectiveValue = 0;
        double EUD0, aux_gEUD, v;
        for (int i = 1; i < o.length; i++) {
            EUD0 = o[i].desiredDose;
            aux_gEUD = this.gEUD[i];
            v = o[i].v;
            objectiveValue = objectiveValue - Math.log(Math.pow(1 + Math.pow(aux_gEUD / EUD0, v), -1));
        }
        this.singleObjectiveValue = objectiveValue;
        System.out.println("Log function Value = " + this.singleObjectiveValue);
    }
    public void getULogFunctionValue(Organs[] o) throws IOException {
        double objectiveValue = 0;
        double EUD0, aux_gEUD, v;
        for (int i = 0; i < o.length; i++) {
            EUD0 = o[i].desiredDose;
            aux_gEUD = this.gEUD[i];
            v = o[i].v;
            if(o[i].isTarget){
                objectiveValue = objectiveValue - Math.log(Math.pow(1 + Math.pow( EUD0 / aux_gEUD , v), -1));
            }else{
                objectiveValue = objectiveValue - Math.log(Math.pow(1 + Math.pow(aux_gEUD / EUD0, v), -1));
            }
        }
        this.singleObjectiveValue = objectiveValue;
        System.out.println("ULog function Value = " + this.singleObjectiveValue);
    }
    public void getUnconstrainedWeightedSum_gEUD_Value(Organs[] o) throws IOException {
        double objectiveValue = 0;
        double EUD0, aux_gEUD, v;
        for (int i = 1; i < o.length; i++) {
            EUD0 = o[i].desiredDose;
            aux_gEUD = this.gEUD[i];
            v = o[i].v;
            if(o[i].isTarget){
                if(aux_gEUD<=EUD0){
                    objectiveValue=objectiveValue+this.weights[i]*Math.pow((EUD0-aux_gEUD),2);
                }else{
                    objectiveValue=objectiveValue+0;
                }    
            }else{
                if(aux_gEUD>=EUD0){
                    objectiveValue=objectiveValue+this.weights[i]*Math.pow((aux_gEUD-EUD0),2);
                }else{
                    objectiveValue=objectiveValue+0;
                }    
            }
        }
        if(this.gEUD[0]>=74){//parametrizar valor 74 tambien en AMPL_Solver
            objectiveValue=objectiveValue+this.weights[0]*Math.pow((this.gEUD[0]-74),2);
        }else{
            objectiveValue=objectiveValue+0;
        }    
        /*
        for (int i = 1; i < o.length; i++) {
            EUD0 = o[i].desiredDose;
            aux_gEUD = this.gEUD[i];
            v = o[i].v;
            objectiveValue = objectiveValue - Math.log(Math.pow(1 + Math.pow(aux_gEUD / EUD0, v), -1));
        }*/
        this.singleObjectiveValue = objectiveValue;
        System.out.println("Log function Value = " + this.singleObjectiveValue);
    }
    public void getVariableWeightsValue(Organs[] o) throws IOException {
        double objectiveValue = 0;
        double aux_gEUD, aux_weights;
        for (int i = 1; i < o.length; i++) {
            aux_gEUD = this.gEUD[i];
            aux_weights = this.weights[i];
            objectiveValue = objectiveValue + (aux_weights * aux_gEUD);
        }
        this.singleObjectiveValue = objectiveValue;
        System.out.println("Weighted Sum with variable weights = " + this.singleObjectiveValue);
    }

    public void getAdaptiveWeightedSumValue(Organs[] o) throws IOException {
        double objectiveValue = 0;
        double aux_gEUD, aux_weights;
        for (int i = 1; i < o.length; i++) {
            aux_gEUD = this.gEUD[i];
            aux_weights = this.weights[i];
            objectiveValue = objectiveValue + (aux_weights * aux_gEUD);
        }
        this.singleObjectiveValue = objectiveValue;
        System.out.println("Adaptive Weighted Sum with variable weights = " + this.singleObjectiveValue);
    }

    public void getConvexLogFunctionValue(Organs[] o) throws IOException {
        double objectiveValue = 0;
        double EUD0, aux_gEUD, v;
        for (int i = 1; i < o.length; i++) {
            EUD0 = o[i].desiredDose;
            aux_gEUD = this.gEUD[i];
            v = o[i].v;
            objectiveValue = objectiveValue - Math.log(Math.pow(Math.pow(aux_gEUD / EUD0, v), -1));
        }
        this.singleObjectiveValue = objectiveValue;
        System.out.println("Inverse function Value = " + this.singleObjectiveValue);
    }

    public void getWeightedSumValue(Organs[] o) throws IOException {
        double objectiveValue = 0;
        double weight, aux_gEUD, v;
        for (int i = 1; i < o.length; i++) {
            if (!o[i].isTarget) {
                weight = o[i].weight;
                aux_gEUD = this.gEUD[i];
                objectiveValue = objectiveValue + aux_gEUD * weight;
            }
        }
        this.singleObjectiveValue = objectiveValue;
        System.out.println("Weighted Sum Value = " + this.singleObjectiveValue);
    }

    public void getVxDx(String jobThreadID, Organs[] o, int[][] valueVx, int[][] valueDx) throws FileNotFoundException, IOException {
        // load dose vector (dose deposited at each voxel
        for (Organs o1 : o) {
            String dir = jobThreadID + "DVH_" + o1.name + ".txt";
            double[] dvh = new double[o1.totalVoxels];
            System.out.println(dir + "\t" + o1.totalVoxels);
            String[] auxReader = null;
            File f = new File(dir);
            if (f.exists()) {
                try(BufferedReader fileIn = new BufferedReader(new FileReader(f))){
                    String line = "";
                    line = fileIn.readLine(); //avoid first line;
                    line = fileIn.readLine();
                    auxReader = line.split(" ");
                    while (!";".equals(auxReader[0])) {
                        for (int j = 0; j < auxReader.length; j++) {
                            if (!"".equals(auxReader[j])) {
                                int auxIndex = (int) Double.parseDouble(auxReader[j]);
                                while ("".equals(auxReader[j + 1])) {
                                    j++;
                                }
                                dvh[auxIndex - 1] = Double.parseDouble(auxReader[j + 1]);
                                j++;
                            }
                        }
                        line = fileIn.readLine();
                        auxReader = line.split(" ");
                    }
                    fileIn.close();
                }
            } else {
                System.out.println("ERROR: CurrentSol file wasn't generated (" + dir + ")");
            }
            Arrays.sort(dvh);
            //calculate Vx indexes
            if (valueVx[o1.index] != null) {
                this.Vx[o1.index] = new double[valueVx[o1.index].length];
                for (int j = 0; j < valueVx[o1.index].length; j++) {
                    double limit = (valueVx[o1.index][j] * o1.desiredDose) / 100;
                    double count = 0;
                    for (int k = o1.totalVoxels - 1; k >= 0; k--) {
                        if (dvh[k] > limit) {
                            count++;
                        } else {
                            break;
                        }
                    }
                    this.Vx[o1.index][j] = (count / o1.totalVoxels) * 100;
                }
            }
            //calculate Dx indexes
            if (valueDx[o1.index] != null) {
                this.Dx[o1.index] = new double[valueDx[o1.index].length];
                for (int j = 0; j < valueDx[o1.index].length; j++) {
                    int limit = o1.totalVoxels - (int) (valueDx[o1.index][j] * o1.totalVoxels) / 100;
                    this.Dx[o1.index][j] = dvh[limit];
                }
            }
        }
    }
   public void getCuadraticSumValue(String jobThreadID,Organs[] o) throws IOException { //Maicholl 
        //Versi贸n Dosimetrica
        double objectiveValue = 0;
        double aux_objective = 0;
        double EUD0, UB;
        this.quadFunction = new double[o.length]; 
        int i=0;
        for (Organs o1 : o) {
            aux_objective = 0;
            EUD0 = o[i].desiredDose;
            String dir = jobThreadID + "DVH_" + o1.name + ".txt";
            double[] intensity = new double[o1.totalVoxels];
            System.out.println(dir + "\t" + o1.totalVoxels);
            String[] auxReader = null;
            File f = new File(dir);
            if (f.exists()) {
                try(BufferedReader fileIn = new BufferedReader(new FileReader(f))){
                    String line = "";
                    line = fileIn.readLine(); //avoid first line;
                    line = fileIn.readLine();
                    auxReader = line.split(" ");
                    while (!";".equals(auxReader[0])) {
                        for (int j = 0; j < auxReader.length; j++) {
                            if (!"".equals(auxReader[j])) {
                                int auxIndex = (int) Double.parseDouble(auxReader[j]);
                                while ("".equals(auxReader[j + 1])) {
                                    j++;
                                }
                                intensity[auxIndex - 1] = Double.parseDouble(auxReader[j + 1]);
                                if(i == 0){
                                    aux_objective = aux_objective + Math.pow(Math.max((EUD0 - intensity[auxIndex - 1]),0),2);
                                }else{
                                    aux_objective = aux_objective + Math.pow(Math.max((intensity[auxIndex - 1] - EUD0),0),2);
                                }
                                j++;
                            }
                        }
                        line = fileIn.readLine();
                        auxReader = line.split(" ");
                    }
                    fileIn.close();
                }
                this.quadFunction[i] = (1/(double)o1.totalVoxels)*aux_objective;//Maicholl
                objectiveValue = objectiveValue + (1/(double)o1.totalVoxels)*aux_objective;
            } else {
                System.out.println("ERROR: CurrentSol file wasn't generated (" + dir + ")");
                objectiveValue = 100000000; //really big number
                break;
            }
            //this.quadFunction[i] = (1/(double)o1.totalVoxels)*aux_objective;//Maicholl  
            //objectiveValue = objectiveValue + (1/(double)o1.totalVoxels)*aux_objective;
            i++;
        }
        //this.quadFunction[0] = o[0].weight*this.quadFunction[0] + o[1].weight*this.quadFunction[1] + o[2].weight*this.quadFunction[2];
        
        this.singleObjectiveValue = objectiveValue;
        System.out.println("Cuadratic Sum Value = " + this.singleObjectiveValue);
       
    }
   
   
    public void getSlope(String jobID, int lineNum) throws FileNotFoundException, IOException {
        String dir = jobID + "lagrangeMultiplier.txt";
        File f = new File(dir);
        try(BufferedReader brLM = new BufferedReader(new FileReader(f))){
            String line = "";
            String[] auxReader;
            for (int i = 0; i < lineNum; i++) {
                line = brLM.readLine();
            }
            auxReader = line.split("  ");
            this.slope = 1 / Double.parseDouble(auxReader[auxReader.length - 1]);
            brLM.close();
        }
    }

    //public void solveLexiRectum(Organs_bkp[] o, DDM_bkp M) throws IOException{
    public void solveLexiRectum(Organs[] o, String jobThreadID) throws IOException {
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver lexRect = new AMPL_Solver(o, this, 0, jobThreadID);

        if (this.intensity.length > 0) {
            lexRect.generateParametersFile(this.intensity);
        } else {
            lexRect.generateParametersFile();
        }
        lexRect.runLexRectum_Solver(o);
        lexRect.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(lexRect.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
        this.getSlope(jobThreadID, 3);
    }
    public void solveLexiRectum(Organs[] o,String solver,String jobThreadID) throws IOException {
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver lexRect = new AMPL_Solver(o, this, solver, 0, jobThreadID);

        if (this.intensity.length > 0) {
            lexRect.generateParametersFile(this.intensity);
        } else {
            lexRect.generateParametersFile();
        }
        lexRect.runLexRectum_Solver(o);
        lexRect.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(lexRect.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
        this.getSlope(jobThreadID, 3);
    }
    //public void solveLexiBladder(Organs_bkp[] o, DDM_bkp M) throws IOException{

    public void solveLexiBladder(Organs[] o, String jobThreadID) throws IOException {
        //double [ ] epsilon_d = {0.0, 0.0, 0.0};
        AMPL_Solver lexBladder = new AMPL_Solver(o, this, 0, jobThreadID);
        if (this.intensity.length > 0) {
            lexBladder.generateParametersFile(this.intensity);
        } else {
            lexBladder.generateParametersFile();
        }
        lexBladder.runLexBladder_Solver(o);
        lexBladder.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(lexBladder.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
        this.getSlope(jobThreadID, 3);
    }

    public void solve_gEUD_Rectum(Organs[] o, double epsilon, String solver, String jobThreadID) throws IOException {
        AMPL_Solver gEUD_Model = new AMPL_Solver(o,this, solver, 30, epsilon, jobThreadID);
        if (this.intensity.length > 0) {
            gEUD_Model.generateParametersFile(this.intensity);
        } else {
            gEUD_Model.generateParametersFile();
        }
        gEUD_Model.run_gEUD_Rectum_Solver(o);
        gEUD_Model.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(gEUD_Model.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
        this.getSlope(jobThreadID, 4);
    }

    public void solve_gEUD_Bladder(Organs[] o, double epsilon, String solver, String jobThreadID) throws IOException {
        AMPL_Solver gEUD_Model = new AMPL_Solver(o,this, solver, 30, epsilon, jobThreadID);
        if (this.intensity.length > 0) {
            gEUD_Model.generateParametersFile(this.intensity);
        } else {
            gEUD_Model.generateParametersFile();
        }
        gEUD_Model.run_gEUD_Bladder_Solver(o);
        gEUD_Model.getSolution();
        this.intensity = new double[this.beamlets];
        System.arraycopy(gEUD_Model.x, 0, this.intensity, 0, this.beamlets);
        //calculate gEUDs
        System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
        this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
        this.getSlope(jobThreadID, 4);
    }

    public double[] getEUD(Organs[] o, String jobThreadID) throws FileNotFoundException, IOException {
        File[] f = new File[o.length];
        BufferedReader[] fileIn = new BufferedReader[o.length];
        double[] aux_gEUD = new double[this.gEUD.length];
        for (int y = 0; y < o.length; y++) {
            f[y] = new File("./" + jobThreadID + "gEUD_" + o[y].name + ".txt");
            fileIn[y] = new BufferedReader(new FileReader(f[y]));
            String line = fileIn[y].readLine();
            line = fileIn[y].readLine();
            String[] auxReader = null;
            if (line != null) {
                auxReader = line.split(" = ");
                aux_gEUD[y] = Double.parseDouble(auxReader[1]);
                //System.out.println(o[y].name +" : " +  aux_gEUD[y] + " [Gy] ");
            } else {
                System.out.println("error: No gEUD value for: " + o[y].name + "(./" + jobThreadID + "gEUD_" + o[y].name + ".txt)");
            }
            fileIn[y].close();
            if (o[y].isTarget) {
                // get gEUD for OAR-PTV
                File g = new File("./" + jobThreadID + "gEUD_" + o[y].name + "_UB.txt");
                try(BufferedReader gIn = new BufferedReader(new FileReader(g))){
                    line = gIn.readLine();
                    line = gIn.readLine();
                    auxReader = null;
                    if (line != null) {
                        auxReader = line.split(" = ");
                        aux_gEUD[this.gEUD.length - 1] = Double.parseDouble(auxReader[1]);
                        //System.out.println(o[y].name +" : " +  aux_gEUD[y] + " [Gy] ");
                    } else {
                        System.out.println("error: No gEUD value for OAR-PTV (./" + jobThreadID + "gEUD_" + o[y].name + "_UB.txt)");
                    }
                    gIn.close();
                }
            }
        }
        return aux_gEUD;
    }

    public boolean sameBAC(Beam[] bac1) {
        boolean flag = false;
        for (int j = 0; j < bac1.length; j++) {
            flag = false; //Assuming bac1 is not equal to bac2
            for (int k = 0; k < this.selAngles.length; k++) {
                if (bac1[j].index == this.selAngles[k].index) {
                    flag = true;
                    break;
                }
            }
            if (flag == false) {
                break;
            }
        }
        return flag;
    }
    
    public void print() {
    	for(double i:scores) {
    		System.out.print(i+" ");
    	}
    	System.out.println("");
    }

    public void printSol(String dirFile) throws IOException {
        BufferedWriter bwFile = null;
        File solFile = new File(dirFile);
        if (solFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(solFile, true));
        } else {
            
            bwFile = new BufferedWriter(new FileWriter(solFile));
            /**************IMPRIME PRIMERA LINEA DEL ARCHIVO*****************/
            for (int i = 0; i < this.beams; i++) {
                writeLine( "Angles["+i+"] \t", bwFile);
            }
            for (int i = 0; i < this.gEUD.length; i++) {
                writeLine( "gEUD["+i+"] \t", bwFile);
            }
            writeLine( "singleObjectiveValue \t", bwFile);
            for (int j = 0; j < this.scores.length; j++) {//katy
              writeLine( "scores["+j+"] \t", bwFile);
            }
            writeLine("globalScore \t", bwFile);//
            for (int j = 0; j < this.quadFunction.length; j++) {//Maichol
              writeLine( "quadFOs["+j+"] \t", bwFile);
            }
            for (int j = 0; j < this.weights.length; j++) {
                writeLine("weights["+j+"] \t", bwFile);
            }
            for (double[] Vx1 : this.Vx) {
                if (Vx1 != null) {
                    for (int j = 0; j < Vx1.length; j++) {
                        writeLine("Vx1["+j+"] \t", bwFile);
                    }
                }
            }
            for (double[] Dx1 : this.Dx) {
                if (Dx1 != null) {
                    for (int j = 0; j < Dx1.length; j++) {
                        writeLine("Dx1["+j+"] \t", bwFile);
                    }
                }
            }
            if (this.allScores != null || this.allScores.length!=0 ){
                for (int i = 0; i < this.allScores.length;i++) {//katy
                    if(this.allScores[i]==null || this.allScores[i].length==0){
                    }else{
                        for (int j = 0; j < this.allScores[i].length; j++) {
                            writeLine("allScores["+i+"]["+j+"] \t", bwFile);
                        } 
                    }
                }
            }
            if (this.Dxx != null || this.Dxx.length!=0 ){
                for (int i = 0; i < this.Dxx.length;i++) {//katy
                    if(this.Dxx[i]==null || this.Dxx[i].length==0){
                    }else{
                        for (int j = 0; j < this.Dxx[i].length; j++) {
                            writeLine("Dxx y Vxx ["+i+"]["+j+"] ", bwFile);
                            switch(j){
                                case 0:
                                    if(i==0){
                                        writeLine("D95", bwFile);
                                    }else{
                                        writeLine("Dmean", bwFile);
                                    }
                                    break;
                                case 1:
                                    if(i==0){
                                        writeLine("D100", bwFile);
                                    }else{
                                        writeLine("Dmax", bwFile);
                                    }
                                    break;
                                case 2:
                                    if(i==0){
                                        writeLine("V95", bwFile);
                                    }else{
                                        writeLine("D50", bwFile);
                                    }
                                    break;
                                case 3:
                                    if(i==0){
                                        writeLine("gEUDbased", bwFile);
                                    }else{
                                        writeLine("D75", bwFile);
                                    }
                                    break;
                                case 4:
                                    if(i==0){
                                        writeLine("H_index", bwFile);
                                    }else{
                                        writeLine("D100", bwFile);
                                    }
                                    break;
                                case 5:
                                    if(i==0){
                                        writeLine("D5", bwFile);
                                    }else{
                                        writeLine("Dmax/Dc", bwFile);
                                    }
                                    break;
                                case 6:
                                    if(i==0){
                                        //writeLine("H_index", bwFile);
                                    }else{
                                        writeLine("Dmean/Dc", bwFile);
                                    }
                                    break;
                            }
                            writeLine("\t", bwFile);
                        } 
                    }
                }
            }
            writeLine( "\n", bwFile);
            /**************IMPRIME PRIMERA LINEA DEL ARCHIVO*****************/
        }
        for (int i = 0; i < this.beams; i++) {
            writeLine(this.selAngles[i].index + "\t", bwFile);
        }
        for (int j = 0; j < this.gEUD.length; j++) {
            writeLine(this.gEUD[j] + "\t", bwFile);
        }
        writeLine(this.singleObjectiveValue + "\t", bwFile);
        
        //katty colocar scores y globalscore y antes de los pesos
        for (int j = 0; j < this.scores.length; j++) {//katy
            writeLine(this.scores[j] + "\t", bwFile);
        }
        writeLine(this.globalScore + "\t", bwFile);//katy
        
        for (int j = 0; j < this.quadFunction.length; j++) {//Maicholl
            writeLine(this.quadFunction[j] + "\t", bwFile);
        }
        
        for (int j = 0; j < this.weights.length; j++) {
            writeLine(this.weights[j] + "\t", bwFile);
        }

        for (double[] Vx1 : this.Vx) {
            if (Vx1 != null) {
                for (int j = 0; j < Vx1.length; j++) {
                    writeLine(Vx1[j] + "\t", bwFile);
                }
            }
        }
        for (double[] Dx1 : this.Dx) {
            if (Dx1 != null) {
                for (int j = 0; j < Dx1.length; j++) {
                    writeLine(Dx1[j] + "\t", bwFile);
                }
            }
        }
        if (this.allScores != null || this.allScores.length!=0 ){
            for (int i = 0; i < this.allScores.length;i++) {//katy
                if(this.allScores[i]==null || this.allScores[i].length==0){
                }else{
                    for (int j = 0; j < this.allScores[i].length; j++) {
                        writeLine(this.allScores[i][j] + "\t", bwFile);
                    }  
                }
            }  
        }
        if (this.Dxx != null || this.Dxx.length!=0 ){
            for (int i = 0; i < this.Dxx.length;i++) {//katy
                if(this.Dxx[i]==null || this.Dxx[i].length==0){
                }else{
                    for (int j = 0; j < this.Dxx[i].length; j++) {
                        writeLine(this.Dxx[i][j] + "\t", bwFile);
                    }  
                }
            }  
        }
        writeLine("x:\t", bwFile);

        for (int j = 0; j < this.intensity.length; j++) {
            writeLine(this.intensity[j] + "\t", bwFile);
        }
        writeLine("\n", bwFile);
        bwFile.close();
    }

    public void printSol(String dirFile, int option) throws IOException {
        BufferedWriter bwFile = null;
        File solFile = new File(dirFile);
        if (solFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(solFile, true));
        } else {
            bwFile = new BufferedWriter(new FileWriter(solFile));
            /**************IMPRIME PRIMERA LINEA DEL ARCHIVO*****************/
            for (int i = 0; i < this.beams; i++) {
                writeLine( "Angles["+i+"] \t", bwFile);
            }
            for (int i = 0; i < this.gEUD.length; i++) {
                writeLine( "gEUD["+i+"] \t", bwFile);
            }
            writeLine( "singleObjectiveValue \t", bwFile);
            for (int j = 0; j < this.scores.length; j++) {//katy
              writeLine( "scores["+j+"] \t", bwFile);
            }
            writeLine("globalScore \t", bwFile);//
            for (int j = 0; j < this.quadFunction.length; j++) {//Maichol
              writeLine( "quadFOs["+j+"] \t", bwFile);
            }
            for (int j = 0; j < this.weights.length; j++) {
                writeLine("weights["+j+"] \t", bwFile);
            }
            if (this.Dxx != null || this.Dxx.length!=0 ){
                for (int i = 0; i < this.Dxx.length;i++) {//katy
                    if(this.Dxx[i]==null || this.Dxx[i].length==0){
                    }else{
                        for (int j = 0; j < this.Dxx[i].length; j++) {
                            writeLine("Dxx y Vxx ["+i+"]["+j+"] ", bwFile);
                            switch(j){
                                case 0:
                                    if(i==0){
                                        writeLine("D95", bwFile);
                                    }else{
                                        writeLine("Dmean", bwFile);
                                    }
                                    break;
                                case 1:
                                    if(i==0){
                                        writeLine("D100", bwFile);
                                    }else{
                                        writeLine("Dmax", bwFile);
                                    }
                                    break;
                                case 2:
                                    if(i==0){
                                        writeLine("V95", bwFile);
                                    }else{
                                        writeLine("D50", bwFile);
                                    }
                                    break;
                                case 3:
                                    if(i==0){
                                        writeLine("gEUDbased", bwFile);
                                    }else{
                                        writeLine("D75", bwFile);
                                    }
                                    break;
                                case 4:
                                    if(i==0){
                                        writeLine("H_index", bwFile);
                                    }else{
                                        writeLine("D100", bwFile);
                                    }
                                    break;
                                case 5:
                                    if(i==0){
                                        writeLine("D5", bwFile);
                                    }else{
                                        writeLine("Dmax/Dc", bwFile);
                                    }
                                    break;
                                case 6:
                                    if(i==0){
                                        //writeLine("H_index", bwFile);
                                    }else{
                                        writeLine("Dmean/Dc", bwFile);
                                    }
                                    break;
                            }
                            writeLine("\t", bwFile);
                        } 
                    }
                }
            }
            for (double[] Vx1 : this.Vx) {
                if (Vx1 != null) {
                    for (int j = 0; j < Vx1.length; j++) {
                        writeLine("Vx1["+j+"] \t", bwFile);
                    }
                }
            }
            for (double[] Dx1 : this.Dx) {
                if (Dx1 != null) {
                    for (int j = 0; j < Dx1.length; j++) {
                        writeLine("Dx1["+j+"] \t", bwFile);
                    }
                }
            }
            if (this.allScores != null || this.allScores.length!=0 ){
                for (int i = 0; i < this.allScores.length;i++) {//katy
                    if(this.allScores[i]==null || this.allScores[i].length==0){
                    }else{
                        for (int j = 0; j < this.allScores[i].length; j++) {
                            writeLine("allScores["+i+"]["+j+"] \t", bwFile);
                        } 
                    }
                }
            }
            writeLine( "\n", bwFile);
            /**************IMPRIME PRIMERA LINEA DEL ARCHIVO*****************/
        }
        for (int i = 0; i < this.beams; i++) {
            writeLine(this.selAngles[i].index + "\t", bwFile);
        }
        for (int j = 0; j < this.gEUD.length; j++) {
            writeLine(this.gEUD[j] + "\t", bwFile);
        }
        writeLine(this.singleObjectiveValue + "\t", bwFile);
        //katty colocar scores y globalscore y antes de los pesos
        for (int j = 0; j < this.scores.length; j++) {//katy
            writeLine(this.scores[j] + "\t", bwFile);
        }
        writeLine(this.globalScore + "\t", bwFile);//katy
        for (int j = 0; j < this.quadFunction.length; j++) {//Maicholl
            writeLine(this.quadFunction[j] + "\t", bwFile);
        }
        for (int j = 0; j < this.weights.length; j++) {
            writeLine(this.weights[j] + "\t", bwFile);
        }

        for (double[] Vx1 : this.Vx) {
            if (Vx1 != null) {
                for (int j = 0; j < Vx1.length; j++) {
                    writeLine(Vx1[j] + "\t", bwFile);
                }
            }
        }
        for (double[] Dx1 : this.Dx) {
            if (Dx1 != null) {
                for (int j = 0; j < Dx1.length; j++) {
                    writeLine(Dx1[j] + "\t", bwFile);
                }
            }
        }
        if (this.allScores != null || this.allScores.length!=0 ){
            for (int i = 0; i < this.allScores.length;i++) {//katy
                if(this.allScores[i]==null || this.allScores[i].length==0){
                }else{
                    for (int j = 0; j < this.allScores[i].length; j++) {
                        writeLine(this.allScores[i][j] + "\t", bwFile);
                    } 
                }
            }
        }
        if (this.Dxx != null || this.Dxx.length!=0 ){
            for (int i = 0; i < this.Dxx.length;i++) {//katy
                if(this.Dxx[i]==null || this.Dxx[i].length==0){
                }else{
                    for (int j = 0; j < this.Dxx[i].length; j++) {
                        writeLine(this.Dxx[i][j] + "\t", bwFile);
                    }  
                }
            }  
        }
        writeLine("\n", bwFile);
        bwFile.close();
    }

    public void printSol(String dirFile, long time) throws IOException {
        BufferedWriter bwFile = null;
        File solFile = new File(dirFile);
        if (solFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(solFile, true));
        } else {
            bwFile = new BufferedWriter(new FileWriter(solFile));
            /**************IMPRIME PRIMERA LINEA DEL ARCHIVO*****************/
            for (int i = 0; i < this.beams; i++) {
                writeLine( "Angles["+i+"] \t", bwFile);
            }
            for (int i = 0; i < this.gEUD.length; i++) {
                writeLine( "gEUD["+i+"] \t", bwFile);
            }
            writeLine( "singleObjectiveValue \t", bwFile);
            writeLine( "time \t", bwFile);
            for (int j = 0; j < this.scores.length; j++) {//katy
              writeLine( "scores["+j+"] \t", bwFile);
            }
            writeLine("globalScore \t", bwFile);//
            for (int j = 0; j < this.weights.length; j++) {
                writeLine("weights["+j+"] \t", bwFile);
            }
            for (double[] Vx1 : this.Vx) {
                if (Vx1 != null) {
                    for (int j = 0; j < Vx1.length; j++) {
                        writeLine("Vx1["+j+"] \t", bwFile);
                    }
                }
            }
            for (double[] Dx1 : this.Dx) {
                if (Dx1 != null) {
                    for (int j = 0; j < Dx1.length; j++) {
                        writeLine("Dx1["+j+"] \t", bwFile);
                    }
                }
            }
            if (this.allScores != null || this.allScores.length!=0 ){
                for (int i = 0; i < this.allScores.length;i++) {//katy
                    if(this.allScores[i]==null || this.allScores[i].length==0){
                    }else{
                        for (int j = 0; j < this.allScores[i].length; j++) {
                            writeLine("allScores["+i+"]["+j+"] \t", bwFile);
                        } 
                    }
                }
            }
            if (this.Dxx != null || this.Dxx.length!=0 ){
                for (int i = 0; i < this.Dxx.length;i++) {//katy
                    if(this.Dxx[i]==null || this.Dxx[i].length==0){
                    }else{
                        for (int j = 0; j < this.Dxx[i].length; j++) {
                            writeLine("Dxx y Vxx ["+i+"]["+j+"] ", bwFile);
                            switch(j){
                                case 0:
                                    if(i==0){
                                        writeLine("D95", bwFile);
                                    }else{
                                        writeLine("Dmean", bwFile);
                                    }
                                    break;
                                case 1:
                                    if(i==0){
                                        writeLine("D100", bwFile);
                                    }else{
                                        writeLine("Dmax", bwFile);
                                    }
                                    break;
                                case 2:
                                    if(i==0){
                                        writeLine("V95", bwFile);
                                    }else{
                                        writeLine("D50", bwFile);
                                    }
                                    break;
                                case 3:
                                    if(i==0){
                                        writeLine("gEUDbased", bwFile);
                                    }else{
                                        writeLine("D75", bwFile);
                                    }
                                    break;
                                case 4:
                                    if(i==0){
                                        writeLine("H_index", bwFile);
                                    }else{
                                        writeLine("D100", bwFile);
                                    }
                                    break;
                                case 5:
                                    if(i==0){
                                        writeLine("D5", bwFile);
                                    }else{
                                        writeLine("Dmax/Dc", bwFile);
                                    }
                                    break;
                                case 6:
                                    if(i==0){
                                        //writeLine("H_index", bwFile);
                                    }else{
                                        writeLine("Dmean/Dc", bwFile);
                                    }
                                    break;
                            }
                            writeLine("\t", bwFile);
                        } 
                    }
                }
            }            
            writeLine( "\n", bwFile);
            /**************IMPRIME PRIMERA LINEA DEL ARCHIVO*****************/
        }
        for (int i = 0; i < this.beams; i++) {
            writeLine(this.selAngles[i].index + "\t", bwFile);
        }
        for (int j = 0; j < this.gEUD.length; j++) {
            writeLine(this.gEUD[j] + "\t", bwFile);
        }
        writeLine(this.singleObjectiveValue + "\t", bwFile);
        writeLine(time + "\t", bwFile);
        //writeLine("x[]:\t", bwFile);
        //katty colocar scores y globalscore y antes de los pesos
        for (int j = 0; j < this.scores.length; j++) {//katy
            writeLine(this.scores[j] + "\t", bwFile);
        }
        writeLine(this.globalScore + "\t", bwFile);//katy
        for (int j = 0; j < this.weights.length; j++) {
            writeLine(this.weights[j] + "\t", bwFile);
        }

        for (double[] Vx1 : this.Vx) {
            if (Vx1 != null) {
                for (int j = 0; j < Vx1.length; j++) {
                    writeLine(Vx1[j] + "\t", bwFile);
                }
            }
        }
        for (double[] Dx1 : this.Dx) {
            if (Dx1 != null) {
                for (int j = 0; j < Dx1.length; j++) {
                    writeLine(Dx1[j] + "\t", bwFile);
                }
            }
        }
        if (this.allScores != null || this.allScores.length!=0 ){
            for (int i = 0; i < this.allScores.length;i++) {//katy
                if(this.allScores[i]==null || this.allScores[i].length==0){
                }else{
                    for (int j = 0; j < this.allScores[i].length; j++) {
                        writeLine(this.allScores[i][j] + "\t", bwFile);
                    } 
                }
            }
        }
        if (this.Dxx != null || this.Dxx.length!=0 ){
            for (int i = 0; i < this.Dxx.length;i++) {//katy
                if(this.Dxx[i]==null || this.Dxx[i].length==0){
                }else{
                    for (int j = 0; j < this.Dxx[i].length; j++) {
                        writeLine(this.Dxx[i][j] + "\t", bwFile);
                    }  
                }
            }  
        }
        //writeLine("x[]:\t", bwFile);
        for (int j = 0; j < this.intensity.length; j++) {
            writeLine(this.intensity[j] + "\t", bwFile);
        }
        writeLine("\n", bwFile);
        bwFile.close();
    }

    public boolean isDominated(TreatmentPlan s) {
        boolean isDominated = true; //means, s dominates 'this'
        //We assume that first index in gEDU  corresponds to the target
        //which is equal for all the solutions
        for (int i = 1; i < this.gEUD.length - 1; i++) {
            if (this.gEUD[i] < s.gEUD[i] * (1 - s.epsilonDominance[i]) ) {
                isDominated = false;
                break;
            }
        }
        return (isDominated);
    }
    public boolean isDominatedEvaluetion(TreatmentPlan s) {
        boolean isDominated = true; //means, s dominates 'this'
        //We assume that first index in gEDU  corresponds to the target
        //which is equal for all the solutions
        for (int i = 0; i < this.scores.length; i++) {
            if (this.scores[i]  < s.scores[i] ) {
                isDominated = false;
                break;
            }
        }
        return (isDominated);
    }
    
    public boolean isDominated(TreatmentPlan s,int optValue) {
        boolean isDominated = true; //means, s dominates 'this'
        //We assume that first index in gEDU  corresponds to the target
        //which is equal for all the solutions
        switch(optValue){
            case 1:
                for (int i = 1; i < this.gEUD.length - 1; i++) {
                    if (this.gEUD[i] - s.epsilonDominance[i] < s.gEUD[i] ) {
                        isDominated = false;
                        break;
                    }
                }
                break;
            case 2:
                /* Aqui va a depender de la funcion objetivo, por ejemplo: Si es la logFunc esta bien que parta desde 1, pero si es la uLogFunc deber铆a partir desde 0*/
                //Asumimos que TODOS los scores se buscan minimizar. Por ejemplo, si usamos TCP, el score deber铆a ser (1-TCP), etc
                for (int i = 0; i < this.scores.length; i++) {
                    if (this.scores[i] - s.epsilonDominance[i] < s.scores[i] ) {
                        isDominated = false;
                        break;
                    }
                }
                break;
            case 9: // Maicholl
                //Asumimos que TODAS las funciones se buscan minimizar. 
                for (int i = 0; i < this.quadFunction.length; i++) {
                    if (this.quadFunction[i] < s.quadFunction[i] * (1+s.epsilonDominance[i])) {
                        isDominated = false;
                        break;
                    }
                }
                break;
            default:
                for (int i = 1; i < this.gEUD.length - 1; i++) {
                    if (this.gEUD[i] - s.epsilonDominance[i]< s.gEUD[i] ) {
                        isDominated = false;
                        break;
                    }
                }
                break;
        }
        
        return (isDominated);
    }
    
    public boolean isDominated(TreatmentPlan s,int optValue, int typeDominance) {
        boolean isDominated = true; //means, s dominates 'this'
        //We assume that first index in gEDU  corresponds to the target
        //which is equal for all the solutions
        switch(optValue){
            case 1:
                for (int i = 1; i < this.gEUD.length - 1; i++) {
                    if (this.gEUD[i] - s.epsilonDominance[i] < s.gEUD[i] ) {
                        isDominated = false;
                        break;
                    }
                }
                break;
            case 2:
                /* Aqui va a depender de la funcion objetivo, por ejemplo: Si es la logFunc esta bien que parta desde 1, pero si es la uLogFunc deber铆a partir desde 0*/
                //Asumimos que TODOS los scores se buscan minimizar. Por ejemplo, si usamos TCP, el score deber铆a ser (1-TCP), etc
                for (int i = 0; i < this.scores.length; i++) {
                    if (this.scores[i] - s.epsilonDominance[i] < s.scores[i] ) {
                        isDominated = false;
                        break;
                    }
                }
                break;
            case 9: // Maicholl
                //Asumimos que TODAS las funciones se buscan minimizar. 
                for (int i = 0; i < this.quadFunction.length; i++) {
                    if(typeDominance == 0){
                        //if (this.quadFunction[i] * (1-this.epsilonDominance[i]) < s.quadFunction[i] ) { // Maicholl porcentual
                        if (this.quadFunction[i] - this.epsilonDominance[i] < s.quadFunction[i] ) {
                            isDominated = false;
                            break;
                        }
                    }else{
                        //if (this.quadFunction[i]  < s.quadFunction[i] * (1-s.epsilonDominance[i])) { // Maicholl porcentual
                        if (this.quadFunction[i]  < s.quadFunction[i] - s.epsilonDominance[i]) {
                            isDominated = false;
                            break;
                        }
                    }
                }
                break;
            default:
                for (int i = 1; i < this.gEUD.length - 1; i++) {
                    if (this.gEUD[i] - s.epsilonDominance[i]< s.gEUD[i] ) {
                        isDominated = false;
                        break;
                    }
                }
                break;
        }
        
        return (isDominated);
    }

    public void computeNewWeights(TreatmentPlan sol, ArrayList<TreatmentPlan> convexHull, double[] bounds, Organs[] o) {
        double[] w = new double[sol.weights.length];
        double[] vertex1 = new double[3];
        double[] vertex2 = new double[3];
        double[] currentPoint = new double[3];
        double[] originalWeights1 = new double[3];
        double[] originalWeights2 = new double[3];
        double distance = 10000, auxDist; //big number
        boolean weightsDone = false;
        Random r = new Random();

        switch (convexHull.size()) {
            case 1://We have only one non-dominated point --> We return same weights
                w[0] = convexHull.get(0).weights[0]; //we don't use this anyway  
                w[1] = convexHull.get(0).weights[1];
                w[2] = convexHull.get(0).weights[2];
                weightsDone = true;
                System.out.println("There is only one solution in the convex hull.");
                System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                        + "New Weights: [" + w[1] + ", " + w[2] + "] ");
                break;
            case 2://We have two points in the convex hull (we might have more non-dominated points though). 
                System.arraycopy(convexHull.get(0).gEUD, 0, vertex1, 0, vertex1.length);
                System.arraycopy(convexHull.get(1).gEUD, 0, vertex2, 0, vertex2.length);
                System.arraycopy(convexHull.get(0).weights, 0, originalWeights1, 0, originalWeights1.length);
                System.arraycopy(convexHull.get(1).weights, 0, originalWeights2, 0, originalWeights2.length);

                //determine whether this solution is supported or not (i.e. it belongs to the convex hull)
                if (sol.sameBAC(convexHull.get(0).selAngles)) {
                    double[] wLeg1 = new double[3];
                    double[] wLeg2 = new double[3];
                    //In this case, wLeg1 is a vertical line, i.e. [1,0]
                    wLeg1[0] = sol.weights[0]; //we don't use this anyway  
                    wLeg1[1] = bounds[1];
                    wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                    wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000000) / 1000000);
                    wLeg1[2] = (double) 1 - wLeg1[1];
                    wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000) / 1000);
                    if (wLeg1[1] + wLeg1[2] != 1) {
                        double dif = 1 - (wLeg1[1] + wLeg1[2]);
                        wLeg1[2] = wLeg1[2] + dif;
                    }
                    //wLeg1[1]= 1; wLeg1[2]= 0;  
                    if (wLeg1[1] + wLeg1[2] != 1) {
                        double dif = 1 - (wLeg1[1] + wLeg1[2]);
                        wLeg1[2] = wLeg1[2] + dif;
                    }

                    wLeg2[0] = sol.weights[0]; //we don't use this anyway
                    //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                    wLeg2[1] = Math.min(bounds[1], ((vertex2[2] - vertex1[2]) / ((vertex1[1] - vertex2[1]) - (vertex1[2] - vertex2[2]))));
                    wLeg2[1] = Math.max(1 - bounds[2], wLeg2[1]);
                    wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                    wLeg2[2] = (double) 1 - wLeg2[2];
                    wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                    if (wLeg2[1] + wLeg2[2] != 1) {
                        double dif = 1 - (wLeg2[1] + wLeg2[2]);
                        wLeg2[2] = wLeg2[2] + dif;
                    }

                    //Since this solution is a supported one (i.e. it belongs to the convex hull), we calculate the
                    //average of the weights of legs 1 and 2
                    w[0] = sol.weights[0]; //we don't use this anyway
                    w[1] = Math.min(bounds[1], (wLeg1[1] + wLeg2[1]) / 2);
                    w[1] = Math.max(1 - bounds[2], w[1]);
                    w[1] = (double) (Math.floor(w[1] * 1000000) / 1000000);
                    w[2] = (double) 1 - w[1];
                    w[2] = (double) (Math.floor(w[2] * 1000000) / 1000000);
                    if (w[1] + w[2] != 1) {
                        double dif = 1 - (w[1] + w[2]);
                        w[2] = w[2] + dif;
                    }
                    System.out.println("This is a Supported Solution (index = 0)");
                    System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "wLegs (" + wLeg1[1] + ", " + wLeg1[2] + ") and "
                            + "(" + wLeg2[1] + ", " + wLeg2[2] + ") "
                            + "New Weights: [" + w[1] + ", " + w[2] + "] ");
                    weightsDone = true;
                } else if (sol.sameBAC(convexHull.get(1).selAngles)) {
                    double[] wLeg1 = new double[3];
                    double[] wLeg2 = new double[3];
                    wLeg2[0] = sol.weights[0]; //we don't use this anyway
                    //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                    wLeg2[2] = bounds[2];
                    wLeg2[2] = Math.max(1 - bounds[1], wLeg2[2]);
                    wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                    wLeg2[1] = (double) 1 - wLeg2[2];
                    wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                    if (wLeg2[1] + wLeg2[2] != 1) {
                        double dif = 1 - (wLeg2[1] + wLeg2[2]);
                        wLeg2[2] = wLeg2[2] + dif;
                    }

                    wLeg1[0] = sol.weights[0]; //we don't use this anyway
                    //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                    wLeg1[1] = Math.min(bounds[1], ((vertex2[2] - vertex1[2]) / ((vertex1[1] - vertex2[1]) - (vertex1[2] - vertex2[2]))));
                    wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                    wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000) / 1000);
                    wLeg1[2] = (double) 1 - wLeg1[2];
                    wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000) / 1000);
                    if (wLeg1[1] + wLeg1[2] != 1) {
                        double dif = 1 - (wLeg1[1] + wLeg1[2]);
                        wLeg1[2] = wLeg1[2] + dif;
                    }

                    //Since this solution is a supported one (i.e. it belongs to the convex hull), we calculate the
                    //average of the weights of legs 1 and 2
                    w[0] = sol.weights[0]; //we don't use this anyway
                    w[1] = Math.min(bounds[1], (wLeg1[1] + wLeg2[1]) / 2);
                    w[1] = Math.max(1 - bounds[2], w[1]);
                    w[1] = (double) (Math.floor(w[1] * 1000) / 1000);
                    w[2] = (double) 1 - w[1];
                    w[2] = (double) (Math.floor(w[2] * 1000) / 1000);
                    if (w[1] + w[2] != 1) {
                        double dif = 1 - (w[1] + w[2]);
                        w[2] = w[2] + dif;
                    }
                    weightsDone = true;
                    System.out.println("This is a Supported Solution (index = 1)");
                    System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "wLegs (" + wLeg1[1] + ", " + wLeg1[2] + ") and "
                            + "(" + wLeg2[1] + ", " + wLeg2[2] + ") "
                            + "New Weights: [" + w[1] + ", " + w[2] + "] ");
                } else { //case for unsupported solutions
                    //Since there is only one leg ('cause there are only two supported solutions) all the
                    //unsupported solutions will use the same weights
                    w[0] = sol.weights[0]; //we don't use this anyway
                    w[1] = Math.min(bounds[1], ((vertex2[2] - vertex1[2]) / ((vertex1[1] - vertex2[1]) - (vertex1[2] - vertex2[2]))));
                    w[1] = Math.max(1 - bounds[2], w[1]);
                    w[1] = (double) (Math.floor(w[1] * 1000000) / 1000000);
                    w[2] = (double) 1 - w[1];
                    w[2] = (double) (Math.floor(w[2] * 1000000) / 1000000);
                    if (w[1] + w[2] != 1) {
                        double dif = 1 - (w[1] + w[2]);
                        w[2] = w[2] + dif;
                    }
                    System.out.println("This is an unsupported Solution, and there is only one segment in the convex hull.");
                    System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "Vertex Weights (" + originalWeights1[1] + ", " + originalWeights1[2] + ") and "
                            + "(" + originalWeights1[1] + ", " + originalWeights1[2] + ") "
                            + "New Weights (segment weights): [" + w[1] + ", " + w[2] + "] ");
                    weightsDone = true;

                }
                break;
            default:
                //First, check whether sol belongs to the convex hull
                for (int i = 0; i < convexHull.size(); i++) {
                    if (sol.sameBAC(convexHull.get(i).selAngles)) {
                        double[] wLeg1 = new double[3];
                        double[] wLeg2 = new double[3];
                        if (i == 0) { //if sol is the first point in the convex hull
                            System.arraycopy(convexHull.get(0).gEUD, 0, vertex1, 0, vertex1.length);
                            System.arraycopy(convexHull.get(1).gEUD, 0, vertex2, 0, vertex2.length);
                            //System.arraycopy(convexHull.get(0).weights, 0, originalWeights1, 0, originalWeights1.length);
                            //System.arraycopy(convexHull.get(1).weights, 0, originalWeights2, 0, originalWeights2.length);

                            wLeg1[0] = sol.weights[0]; //we don't use this anyway
                            wLeg1[1] = bounds[1];
                            wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                            wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000000) / 1000000);
                            wLeg1[2] = (double) 1 - wLeg1[1];
                            wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000000) / 1000000);
                            if (wLeg1[1] + wLeg1[2] != 1) {
                                double dif = 1 - (wLeg1[1] + wLeg1[2]);
                                wLeg1[2] = wLeg1[2] + dif;
                            }

                            wLeg2[0] = sol.weights[0]; //we don't use this anyway
                            wLeg2[1] = Math.min(bounds[1], ((vertex2[2] - vertex1[2]) / ((vertex1[1] - vertex2[1]) - (vertex1[2] - vertex2[2]))));
                            wLeg2[1] = Math.max(1 - bounds[2], wLeg2[1]);
                            wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                            wLeg2[2] = (double) 1 - wLeg2[2];
                            wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                            if (wLeg2[1] + wLeg2[2] != 1) {
                                double dif = 1 - (wLeg2[1] + wLeg2[2]);
                                wLeg2[2] = wLeg2[2] + dif;
                            }
                            System.out.println("This is a Supported Solution (index = 0)");

                        } else if (i == convexHull.size() - 1) {//if sol is the last point in the convex hull    
                            System.arraycopy(convexHull.get(convexHull.size() - 2).gEUD, 0, vertex1, 0, vertex1.length);
                            System.arraycopy(convexHull.get(convexHull.size() - 1).gEUD, 0, vertex2, 0, vertex2.length);
                            //System.arraycopy(convexHull.get(0).weights, 0, originalWeights1, 0, originalWeights1.length);

                            wLeg2[0] = sol.weights[0]; //we don't use this anyway
                            //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                            wLeg2[2] = bounds[2];
                            wLeg2[2] = Math.max(1 - bounds[1], wLeg2[2]);
                            wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                            wLeg2[1] = (double) 1 - wLeg2[2];
                            wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                            if (wLeg2[1] + wLeg2[2] != 1) {
                                double dif = 1 - (wLeg2[1] + wLeg2[2]);
                                wLeg2[2] = wLeg2[2] + dif;
                            }

                            wLeg1[0] = sol.weights[0]; //we don't use this anyway
                            //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                            wLeg1[1] = Math.min(bounds[1], ((vertex2[2] - vertex1[2]) / ((vertex1[1] - vertex2[1]) - (vertex1[2] - vertex2[2]))));
                            wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                            wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000000) / 1000000);
                            wLeg1[2] = (double) 1 - wLeg1[2];
                            wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000000) / 1000000);
                            if (wLeg1[1] + wLeg1[2] != 1) {
                                double dif = 1 - (wLeg1[1] + wLeg1[2]);
                                wLeg1[2] = wLeg1[2] + dif;
                            }
                            System.out.println("This is a Supported Solution (index = " + (convexHull.size() - 1) + ")");
                        } else { //sol belongs to the convex hull but it is neither the first nor the last point in the convex hull

                            System.arraycopy(convexHull.get(i - 1).gEUD, 0, vertex1, 0, vertex1.length);
                            System.arraycopy(convexHull.get(i).gEUD, 0, currentPoint, 0, currentPoint.length);
                            System.arraycopy(convexHull.get(i + 1).gEUD, 0, vertex2, 0, vertex2.length);
                            //System.arraycopy(convexHull.get(i-1).weights, 0, originalWeights1, 0, originalWeights1.length);
                            //System.arraycopy(convexHull.get(i+1).weights, 0, originalWeights2, 0, originalWeights2.length);

                            wLeg1[0] = sol.weights[0]; //we don't use this anyway
                            wLeg1[1] = Math.min(bounds[1], ((currentPoint[2] - vertex1[2]) / ((vertex1[1] - currentPoint[1]) - (vertex1[2] - currentPoint[2]))));
                            wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                            wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000000) / 1000000);
                            wLeg1[2] = (double) 1 - wLeg1[1];
                            wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000000) / 1000000);
                            if (wLeg1[1] + wLeg1[2] != 0) {
                                double dif = 1 - (wLeg1[1] + wLeg1[2]);
                                wLeg1[2] = wLeg1[2] + dif;
                            }

                            wLeg2[0] = sol.weights[0]; //we don't use this anyway
                            wLeg2[1] = Math.min(bounds[1], ((vertex2[2] - currentPoint[2]) / ((currentPoint[1] - vertex2[1]) - (currentPoint[2] - vertex2[2]))));
                            wLeg2[1] = Math.max(1 - bounds[2], wLeg2[1]);
                            wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                            wLeg2[2] = (double) 1 - wLeg2[1];
                            wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                            if (wLeg2[1] + wLeg2[2] != 0) {
                                double dif = 1 - (wLeg2[1] + wLeg2[2]);
                                wLeg2[2] = wLeg2[2] + dif;
                            }
                            System.out.println("This is a Supported Solution (index = " + i + ")");
                        }
                        //Since this solution is a supported one (i.e. it belongs to the convex hull), we calculate the
                        //average of the weights of legs 1 and 2

                        w[0] = sol.weights[0]; //we don't use this anyway
                        w[1] = Math.min(bounds[1], (wLeg1[1] + wLeg2[1]) / 2);
                        w[1] = Math.max(1 - bounds[2], w[1]);
                        w[1] = (double) (Math.floor(w[1] * 1000000) / 1000000);
                        w[2] = (double) 1 - w[1];
                        w[2] = (double) (Math.floor(w[2] * 1000000) / 1000000);
                        if (w[1] + w[2] != 1) {
                            double dif = 1 - (w[1] + w[2]);
                            w[2] = w[2] + dif;
                        }
                        System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                                + "wLegs (" + wLeg1[1] + ", " + wLeg1[2] + ") and "
                                + "(" + wLeg2[1] + ", " + wLeg2[2] + ") "
                                + "New Weights: [" + w[1] + ", " + w[2] + "] ");
                        weightsDone = true;
                    }
                }
                if (!weightsDone) { //means this solution (sol) in unsuported (i.e. it does not belong to the convex hull)
                    //locate the closest segment of the convex hull to sol
                    for (int i = 0; i < convexHull.size() - 1; i++) {
                        double[] p0 = new double[3];
                        double[] p1 = new double[3];
                        double[] p2 = new double[3];
                        System.arraycopy(convexHull.get(i).gEUD, 0, p0, 0, p0.length);
                        System.arraycopy(convexHull.get(i + 1).gEUD, 0, p1, 0, p1.length);
                        System.arraycopy(sol.gEUD, 0, p2, 0, p2.length);
                        auxDist = Math.abs((p2[2] - p1[2]) * p0[1] - (p2[1] - p1[1]) * p0[2] + p2[1] * p1[2] - p2[2] * p1[1]);
                        auxDist = auxDist / Math.sqrt(Math.pow(p2[2] - p1[2], 2) + Math.pow(p2[1] - p1[1], 2));
                        if (auxDist < distance) {
                            distance = auxDist;
                            System.arraycopy(convexHull.get(i).gEUD, 0, vertex1, 0, vertex1.length);
                            System.arraycopy(convexHull.get(i + 1).gEUD, 0, vertex2, 0, vertex1.length);
                            System.arraycopy(convexHull.get(i).weights, 0, originalWeights1, 0, originalWeights1.length);
                            System.arraycopy(convexHull.get(i + 1).weights, 0, originalWeights2, 0, originalWeights2.length);
                        }
                    }
                    w[0] = sol.weights[0]; //we don't use this anyway
                    w[1] = Math.min(bounds[1], ((vertex2[2] - vertex1[2]) / ((vertex1[1] - vertex2[1]) - (vertex1[2] - vertex2[2]))));
                    w[1] = Math.max(1 - bounds[2], w[1]);
                    w[1] = (double) (Math.floor(w[1] * 1000000) / 1000000);
                    w[2] = (double) 1 - w[1];
                    w[2] = (double) (Math.floor(w[2] * 1000000) / 1000000);
                    if (w[1] + w[2] != 0) {
                        double dif = 1 - (w[1] + w[2]);
                        w[2] = w[2] + dif;
                    }
                    System.out.println("This is an unsupported Solution.");
                    System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "Vertex Weights (" + originalWeights1[1] + ", " + originalWeights1[2] + ") and "
                            + "(" + originalWeights1[1] + ", " + originalWeights1[2] + ") "
                            + "New Weights (segment weights): [" + w[1] + ", " + w[2] + "] ");
                }
                break;
        }
        for (int i = 0; i < this.weights.length; i++) {
            this.weights[i] = w[i];
        }
    }

    public static void writeLine(String l, BufferedWriter bw) throws IOException {
        bw.write(l);
    }
    public static Comparator<TreatmentPlan> scoreRectum_Comparator = new Comparator<TreatmentPlan>() {
        public int compare(TreatmentPlan tp1, TreatmentPlan tp2) {
            if (tp1.scores[1] < tp2.scores[1]) {
                return -1;
            } else {
                return 1;
            }
        }
    };
    public static Comparator<TreatmentPlan> gEUDRectum_Comparator = new Comparator<TreatmentPlan>() {
        public int compare(TreatmentPlan tp1, TreatmentPlan tp2) {
            if (tp1.gEUD[1] < tp2.gEUD[1]) {
                return -1;
            } else {
                return 1;
            }
        }
    };

    public static Comparator<TreatmentPlan> judgementFunct_Comparator = new Comparator<TreatmentPlan>() {
        public int compare(TreatmentPlan tp1, TreatmentPlan tp2) {
            if (tp1.singleObjectiveValue < tp2.singleObjectiveValue) {
                return -1;
            } else {
                return 1;
            }
        }
    };
}
