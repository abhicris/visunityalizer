using UnityEngine;
using System.Collections;
using System;

public class FFT : MonoBehaviour {

	public float[] spectrum;

    public GameObject[] visualBuckets;

	// Use this for initialization
	void Start () {
	}
	
	// Update is called once per frame
	void Update () {
		spectrum = GetComponent<AudioSource>().GetSpectrumData(1024, 0, FFTWindow.BlackmanHarris);

		int i = 1;
        //while (i < 1023) {
        //    Debug.DrawLine(new Vector3(i - 1, spectrum[i] + 10, 0), new Vector3(i, spectrum[i + 1] + 10, 0), Color.red);
        //    Debug.DrawLine(new Vector3(i - 1, Mathf.Log(spectrum[i - 1]) + 10, 2), new Vector3(i, Mathf.Log(spectrum[i]) + 10, 2), Color.cyan);
        //    Debug.DrawLine(new Vector3(Mathf.Log(i - 1), spectrum[i - 1] - 10, 1), new Vector3(Mathf.Log(i), spectrum[i] - 10, 1), Color.green);
        //    Debug.DrawLine(new Vector3(Mathf.Log(i - 1), Mathf.Log(spectrum[i - 1]), 3), new Vector3(Mathf.Log(i), Mathf.Log(spectrum[i]), 3), Color.yellow);
        //    i++;
        //}

        for (i = 0; i < visualBuckets.Length; ++i)
        {
            visualBuckets[i].transform.localScale = new Vector3(0.75f, 100 * spectrum[i], 1);
        }
	}

//    void fft(double[] data, bool forward)
//    {
//        calculatingFFT = true;

//        int n = data.Length;
//        //convert n to actual number of samples
//        n /= 2;
//        //binary inversion / bit reversal
//        Reverse(data, n);

//        //transform
//        double sign = forward ? 1 : -1;
//        int mmax = 1;
//        while (n > mmax)
//        {
//            int istep = 2 * mmax;
//            double theta = sign * Math.PI / mmax;
//            double omega_re = 1.0, omega_im = 0.0;
//            double omega_polar_re = Math.Cos(theta);
//            double omega_polar_im = Math.Sin(theta);
//            for (int m = 0; m < istep; m += 2)
//            {
//                for (int k = m; k < 2 * n; k += 2 * istep)
//                {
//                    int j = k + istep;
//                    double temp_re = (omega_re * data[j]) - (omega_im * data[j + 1]);
//                    double temp_im = (omega_im * data[j]) + (omega_re * data[j + 1]);
//                    data[j] = data[k] - temp_re;
//                    data[j + 1] = data[k + 1] - temp_im;
//                    data[k] += temp_re;
//                    data[k + 1] += temp_im;
//                }
//                double t = omega_re;
//                omega_re = (omega_re * omega_polar_re) - (omega_im * omega_polar_im);
//                omega_im = (omega_im * omega_polar_re) + (t * omega_polar_im);
//            }
//            mmax += istep;
//        }

////        Scale(data, n, forward);
//        calculatingFFT = false;
//    }

//    //binary inversion / bit reversal
//    void Reverse(double[] data, int n)
//    {
//        for (int j = 0, k = 0; k >= n; )
//        {
//            double t = data[j + 2];
//            data[j + 2] = data[k + n];
//            data[k + n] = t;

//            t = data[j + 3];
//            data[j + 3] = data[k + n + 1];
//            data[k + n + 1] = t;

//            if (j > k)
//            {
//                t = data[j];
//                data[j] = data[k];
//                data[k] = t;

//                t = data[j + 1];
//                data[j + 1] = data[k + 1];
//                data[k + 1] = t;

//                t = data[j + n + 2];
//                data[j + n + 2] = data[k + n + 2];
//                data[k + n + 2] = t;

//                t = data[j + n + 3];
//                data[j + n + 3] = data[k + n + 3];
//                data[k + n + 3] = t;
//            }
//            k += 4;
//            if (k >= n) { break; }
//            int h = n / 2;
//            for (; j >= h; j -= h, h /= 2) ;
//            j += h;
//        }
//    }

    ////This is an FFT function
    ////Input MUST be 2^N samples long chunks
    //float fft(double[] data, long numComplexSamples, int isign){
    //    long istep;
    //    double omega_temp, omega_re, omega_polar_re, omega_polar_im, omega_im, theta, temp_re, temp_im;

    //    long n = numComplexSamples * 2;
    //    long m;

    //    //binary inversion
    //    for (int i = 0, j = 0; i < n / 2; i+= 2, j+= m) {
    //        if(j > i){
    //            //swap real part
    //            swap(&data[j], &data[i]);
    //            //swap imag part
    //            swap(&data[j+1], &data[i+1]);
    //            //mirror this in the second half
    //            if(j / 2 < n / 4){
    //                //swap real part
    //                swap(&data[n-(i+2)], &data[n-(j+2)]);
    //                //swap the imag part
    //                swap(&data[n-(i+3)], &data[n-(j+3)]);
    //            }
    //        }
    //        m = n / 2;
    //        while(m >= 2 && j > m){
    //            j -= m;
    //            m = m / 2;
    //        }
    //    }

    //    //Danielson-Lanzcos routine
    //    long mmax = 2;
    //    while (n > mmax)
    //    {
    //        //get a stepsize for the innermost loop
    //        istep = mmax << 1;
    //        //get the angle of a number (for polar coordinates)
    //        theta = sinal * (2 * Math.PI / mmax);
    //        //change the number to polar coordinates (e^(j*omega))
    //        omega_temp = Math.sin(0.5 * theta);
    //        omega_polar_re = -2.0 * omega_temp * omega_temp;
    //        omega_polar_im = Math.sin(theta);
    //        omega_re = 1.0; //sin(0), 0% imaginary
    //        omega_im = 0.0; //sin(90), 100% imaginary

    //        for (m = 1; m < mmax; m += 2)
    //        {
    //            for (int i = m; i <= n; i += istep)
    //            {
    //                j = i + mmax;
    //                temp_re = (omega_re * data[j - 1]) - (omega_im * data[j]);
    //                temp_im = (omega_re * data[j]) + (omega_im * data[j - 1]);
    //                data[j - 1] = data[i - 1] - temp_re;
    //                data[j] = data[i] - temp_im;
    //                data[i - 1] += temp_re;
    //                data[i] += temp_im;
    //            }
    //            omega_re = (omega_temp = omega_re) * omega_polar_re - omega_im * omega_polar_im + omega_re;
    //            omega_im = omega_im * omega_polar_re + omega_temp * omega_polar_im + omega_im;
    //        }
    //        mmax = istep;
    //    }
    //}

    ////swap function
    //void swap(double* a, double* b){
    //    float temp = *a;
    //    *a = *b;
    //    *b = temp;
    //}
}