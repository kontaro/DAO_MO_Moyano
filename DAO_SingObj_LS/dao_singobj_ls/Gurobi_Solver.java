package dao_singobj_ls;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import IMRT_Base.DDM;
import IMRT_Base.Organs;
import IMRT_Base.TreatmentPlan;
import gurobi.*;


/**
*
* @author Moyano
*/

public class Gurobi_Solver {
	
    public int organs;//number of organ
    public int beams;//number of organ
    public int []R;//number of voxels by region
    public int aperture;//number of aperture
    
    public int[] bmlts;// #number of beamlets by angles
    public double[] weight; //weights of objective
    public double[] EUD0;
    
    
    public int totalBmlts; //total beamlets
    public int[] angles;
    
    public double[][] newIntensity; 
    
    public double[] LB;
    public double[] UB;
    public Boolean[] isTarget;
    public double epsilon;
    //public double t;
    public double[] x;
    public String jobThreadID;
    public String solver;
    public double maxIntensity;
    public double minIntensity;
    public double objVal;
    
    public int[] eud;
	GRBEnv env;
	GRBModel model;
	DDM M;
	TreatmentPlan sol;
	public Gurobi_Solver(TreatmentPlan sol,DDM M,Organs[] o,int[] selAngles,int[] dd,double[] weight) throws GRBException {
		beams=sol.beams;
		this.eud=dd;
		R=new int[o.length];
		organs=o.length;
		bmlts=new int[beams];
		this.weight=weight;
		bmlts=selAngles;
		aperture=sol.apertures.get(0).size();
		minIntensity=0;
		maxIntensity=sol.max_intesity;
		this.M=M;
		this.sol=sol;
		for(int i=0;i<R.length;i++) {
			R[i]=M.OrganVoxels[i];
		}
		setEnv();
		setModel();
		//writeModel();
		//model.dispose();
        //env.dispose();
		
	}
	
	public void setNewModel(TreatmentPlan sol) {
		aperture=sol.apertures.get(0).size();
		this.sol=sol;
		boolean error=true;
		do {
		
			try {
				setEnv();
				setModel();
				//model.dispose();
		        //env.dispose();
		        error=false;
			} catch (GRBException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}while(error);
		
		
		
	}
	
	public void reset() {
		try {
			model.reset();
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public void setEnv() throws GRBException {
	      this.env = new GRBEnv(true);
	      env.set("OutputFlag", "0");
	      //env.set("logFile", "mip1.log");
	      
	      env.set(GRB.IntParam.Threads, 4);
	      env.set("Aggregate","0");
	      env.set("Method", "2");
	      env.set("NormAdjust", "2");
	      env.start();
	}
	public void setModel() throws GRBException {
		this.model = new GRBModel(env);
		model.set(GRB.StringAttr.ModelName, "Direct Aperture Optimization");
		
		
		//set Variables
		
		// set variables for intensity an V
		GRBVar[][] intensity= new GRBVar[beams][aperture];//intensity for beam,for aperture
		GRBVar[][] voxel=new GRBVar[R.length][];
		int indexi=0,indexj=0;
		
		for (int i = 0; i < this.beams; ++i) {
            for (int j = 0; j < this.aperture; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	intensity[i][j] = model.addVar(minIntensity,maxIntensity,0.0, GRB.CONTINUOUS ,"Intensity" +"["+ indexi + "." + indexj+"]");
            	intensity[i][j].set(GRB.DoubleAttr.Start, sol.apertures.get(i).get(j).intensity); 
            }
        }
		
		
		for (int i = 0; i < this.R.length; ++i) {
			voxel[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	voxel[i][j] = model.addVar(0.0,50.0,0.0, GRB.CONTINUOUS ,"v" + indexi + "[" + indexj+"]");
                
            }
        }
		

		
		//set constraints
		//for 
	
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Integer key, beamblet, totalBeamblets,beamIndex, count_voxel;
		Double radiation, coefficent;
		int diffBeamblets=0;
		
		for (int o = 0;  o< organs; o++) {
			aux_index = index_dao_ddm.get(o);
			aux_values = value_dao_ddm.get(o);
			keys = aux_index.keys();
			
			//Recorremos claves de voxel por organo para su evaluación
			
			count_voxel = 0;
			while(keys.hasMoreElements()){
				GRBLinExpr voxelRadiation= new GRBLinExpr();
				key = keys.nextElement();
	            beams = aux_index.get(key);
	            
	            for(int b=0;b<beams.size(); b++){ // de aqui vamos a sacar el beam (indice del angulo)
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	beamblet = beams.get(b);
	            	totalBeamblets = 0;
	            	beamIndex = 0;
	            	diffBeamblets=0;
	            	for(int z=0; z<bmlts.length;z++) {
	            		totalBeamblets+= bmlts[z];
	            		if(beamblet < totalBeamblets){
	            			beamIndex = z;
	            			break;
	            		}
	            		diffBeamblets+=bmlts[z];
	            	}
	            	for(int a=0;a<aperture;a++) {
	            		coefficent = (double)sol.aperturesBmlts.get(beamIndex)[a][beamblet-diffBeamblets]*radiation;
	            		//coefficent = (double)aper[beamIndex][a][(beamblet-diffBeamblets)]*radiation;
	            		if(coefficent!=0 && o==0) {
							coefficent=coefficent*-1;
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][a] );
							//System.out.println(coefficent);
						}
						else {
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][a] );
						}
	            		
	            		
	            	}
	            }
	            if(o==0) {
					voxelRadiation.addConstant(eud[o]);
				}
				else{
					int constEud=eud[o]*-1;
					voxelRadiation.addConstant(constEud);
				}
				GRBLinExpr V=new GRBLinExpr();
				V.addTerm(1,voxel[o][count_voxel]);
				model.addConstr(V, GRB.GREATER_EQUAL, voxelRadiation, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				count_voxel++;
			}
		}
		
		//set model
		GRBQuadExpr objFunc= new GRBQuadExpr();
		for(int o = 0;  o< organs; o++) {
			double coef=(double)((weight[o]/R[o]));
			
			for (int j = 0; j < R[o]; ++j) {
				objFunc.addTerm(coef, voxel[o][j],voxel[o][j]);
				System.out.print("");
			}
		}
		
		model.setObjective(objFunc, GRB.MINIMIZE);
		model.update();
        model.optimize();
        model.update();
		double[][]getIntensity=new double[this.beams][this.aperture];
        for (int i = 0; i < this.beams; ++i) {
            for (int j = 0; j < this.aperture; ++j) {
            	getIntensity[i][j]=intensity[i][j].get(GRB.DoubleAttr.X);
            	//String varName="intensity"+i+"."+j;
            	//getIntensity[i][j]=model.getVarByName(varName).get(GRB.DoubleAttr.X);
            	//intensity[i][j] = model.addVar(minIntensity,maxIntensity,0.0, GRB.CONTINUOUS ,"Intensity" + i + "." + j);
                
            }
        }
		//GRBVar[] vars = model.getVars();
        newIntensity=getIntensity;
        objVal=model.get(GRB.DoubleAttr.ObjVal);
        
		
		
	}
	public void writeModel() throws GRBException {
		model.write("out1.lp");
		model.write("out1.mst");
		//model.write("out1.mps");
		model.write("out1.sol");
		
	}
	
	
	public void cuadraticSum(TreatmentPlan sol) {
		
	}
	
	public void returnIntensity() {
		
	}
}
