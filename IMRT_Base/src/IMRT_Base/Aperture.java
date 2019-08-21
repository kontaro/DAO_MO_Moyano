package IMRT_Base;


import java.util.Vector;

public class Aperture {
	
	public int[][] aperture;
	public double intensity;
	
	public Aperture(int[][] a, int i){
		this.aperture = a;
		this.intensity = i;
	}
	public Aperture(Aperture original){
		this.aperture = new int[original.aperture.length][original.aperture[0].length];
		this.intensity = original.intensity;
		
		for(int i=0;i<original.aperture.length;i++) {
			for(int j=0;j<original.aperture[0].length;j++) {
				aperture[i][j]= original.aperture[i][j];
			}
		}
		
		
		
	}
	
	/**
	 * funcion que abre o cierra la apertura(no considera si existen movimientos infactibles)
	 * @param row 
	 * @param state si es 1 se abre la apertura, en caso contrario se cierra
	 * @param direction indica en que direccion se hace la operacion, si es true se hace en la izquierda, si es false en la derecha
	 */
	public void changeRow(int row, int state,boolean direction) {
		int aux=-1;
		if(state==1) aux=1;
		if(direction==true) {
			aperture[row][0]=aperture[row][0]+aux;
		}
		else {
			aperture[row][1]=aperture[row][1]+aux;
		}
	}
}
