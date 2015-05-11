using UnityEngine;
using System.Collections;

public class FFT : MonoBehaviour {

	// Use this for initialization
	void Start () {
	
	}
	
	// Update is called once per frame
	void Update () {
	
	}
	
	//This is an FFT function
	float fft(float[] data, long numComplexSamples, int isign){
		long mmax, istep;
		double omega_temp, omega_re, omega_pi_re, omega_pi_im, omega_im, theta, temp_re, temp_im;

		long n = numComplexSamples * 2;
		long m;

		//binary inversion method
		for (int i = 0, j = 0; i < n / 2; i+= 2, j+= m) {
			if(j > i){
				//swap real part
				swap(&data[j], &data[i]);
				//swap imag part
				swap(&data[j+1], &data[i+1]);
				//mirror this in the second half
				if(j / 2 < n / 4){
					//swap real part
					swap(&data[n-(i+2)], &data[n-(j+2)]);
					//swap the imag part
					swap(&data[n-(i+3)], &data[n-(j+3)]);
				}
			}
			m = n / 2;
			while(m >= 2 && j > m){
				j -= m;
				m = m / 2;
			}
		}
	}

	//swap function
	void swap(float* a, float* b){
		float temp = *a;
		*a = *b;
		*b = temp;
	}
}