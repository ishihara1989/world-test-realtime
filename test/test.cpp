//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//
// test.exe input.wav outout.wav f0 spec
// input.wav  : Input file
//
// output.wav : Output file
// f0         : F0 scaling (a positive number)
// spec       : Formant scaling (a positive number)
//
// Note: This version output three speech synthesized by different algorithms.
//       When the filename is "output.wav", "01output.wav", "02output.wav" and
//       "03output.wav" are generated. They are almost all the same.
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <chrono>

#if (defined (__WIN32__) || defined (_WIN32)) && !defined (__MINGW32__)
#include <conio.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
#pragma warning(disable : 4996)
#endif
#if (defined (__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
#include <stdint.h>
#include <sys/time.h>
#endif

// For .wav input/output functions.
#include "audioio.h"

// WORLD core functions.
// Note: win.sln uses an option in Additional Include Directories.
// To compile the program, the option "-I $(SolutionDir)..\src" was set.
#include "world/d4c.h"
#include "world/dio.h"
#include "world/harvest.h"
#include "world/matlabfunctions.h"
#include "world/cheaptrick.h"
#include "world/stonemask.h"
#include "world/synthesis.h"
#include "world/synthesisrealtime.h"


//-----------------------------------------------------------------------------
// struct for WORLD
// This struct is an option.
// Users are NOT forced to use this struct.
//-----------------------------------------------------------------------------
typedef struct {
  double frame_period;
  int fs;

  double *f0;
  double *time_axis;
  int f0_length;

  double **spectrogram;
  double **aperiodicity;
  int fft_size;
} WorldParameters;

using namespace std;

void convert(double f0mul, double spmul, double* f0, double** sp, double** ap, int f0_length, int hfft, double* buffer){
  for(int i=0; i<f0_length; i++){
    f0[i]*=f0mul;
    //sp
    for(int j=0; j<hfft; j++){
      buffer[j]=log(sp[i][j]+1e-30);
    }
    for(int j=0; j<hfft; j++){
      double target_index = j/spmul;
      if(target_index>=hfft) {
        for(int k=j; k<hfft; k++){
          sp[i][j]=1e-30;
        }
        break;
      }
      int z=(int)target_index;
      double r=target_index-z;
      sp[i][j] = exp((1-r)*buffer[z]+r*buffer[z+1]);
    }
    //ap
    for(int j=0; j<hfft; j++){
      buffer[j]=log(ap[i][j]+1e-30);
    }
    for(int j=0; j<hfft; j++){
      double target_index = j/spmul;
      if(target_index>=hfft) {
        for(int k=j; k<hfft; k++){
          ap[i][j]=1.0;
        }
        break;
      }
      int z=(int)target_index;
      double r=target_index-z;
      ap[i][j] = exp((1-r)*buffer[z]+r*buffer[z+1]);
    }
  }
}

struct OnlineWorldParams{
  double* f0;
  double* time_axis;
  double** spectrogram;
  double** aperiodicity;
};

void online(const double *x, int x_length, int fs, double *y){
  double frame_period = 4.0;
  double f0_floor = 80.0;
  double c_f0_floor = 71.0;
  HarvestOption h_option = { 0 };
  HarvestBuffer h_buffer;
  DioOption dio_option = {0};
  CheapTrickOption c_option = {0};
  D4COption d_option = {0};
  InitializeHarvestOption(&h_option);
  h_option.frame_period = frame_period;
  h_option.f0_floor = f0_floor;
  InitializeDioOption(&dio_option);
  dio_option.frame_period = 4.0;
  dio_option.speed = 1;
  dio_option.f0_floor = 40.0;
  dio_option.allowed_range = 0.1;
  InitializeCheapTrickOption(fs, &c_option);
  c_option.f0_floor = c_f0_floor;
  int fft_size = GetFFTSizeForCheapTrick(fs, &c_option);
  c_option.fft_size = fft_size;
  InitializeD4COption(&d_option);
  d_option.threshold = 0.85;

  WorldSynthesizer synthesizer = { 0 };
  int buffer_size = 64;

  int overlap = 2048;
  int shift=2048;
  int frame_shift_samp = static_cast<int>(frame_period*fs/1000);
  int win_length=shift+overlap;
  int f0_length = GetSamplesForHarvest(fs, win_length, h_option.frame_period);
  int hfft = fft_size/2+1;
  int cut = overlap/2/frame_shift_samp;
  int ring_buffer_size = 4*f0_length;
  double *f0 = new double[ring_buffer_size];
  double *pre_f0 = new double[ring_buffer_size];
  double *time_axis = new double[f0_length];
  double **spectrogram = new double *[ring_buffer_size];
  double **aperiodicity = new double *[ring_buffer_size];
  for(int i=0; i<ring_buffer_size; i++){
    spectrogram[i] = new double[hfft];
    aperiodicity[i] = new double[hfft];
  }

  InitializeHarvestBuffer(win_length, fs, &h_option, &h_buffer);

  double* interp_buffer = new double[hfft];
  
  InitializeSynthesizer(fs, frame_period, fft_size, buffer_size, ring_buffer_size, &synthesizer);

  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
  int index = 0;
  int ring_buffer_index = 0;
  for(int i=0; i<x_length-win_length;i+=shift){
    // HarvestNoMemory(&x[i], win_length, fs, &h_option, time_axis, &f0[ring_buffer_index], &h_buffer);
    // Harvest(&x[i], win_length, fs, &h_option, time_axis, &f0[ring_buffer_index]);
    Dio(&x[i], win_length, fs, &dio_option, time_axis, &pre_f0[ring_buffer_index]);
    StoneMask(x, win_length, fs, time_axis,
        &pre_f0[ring_buffer_index], f0_length, &f0[ring_buffer_index]);
    CheapTrick(&x[i], win_length, fs, time_axis, &f0[ring_buffer_index], f0_length, &c_option, &spectrogram[ring_buffer_index]);
    D4C(&x[i], win_length, fs, time_axis, &f0[ring_buffer_index], f0_length, fft_size, &d_option, &aperiodicity[ring_buffer_index]);

    convert(2.0,1.2,&f0[ring_buffer_index],&spectrogram[ring_buffer_index],&aperiodicity[ring_buffer_index],f0_length,hfft,interp_buffer);
    for (int ii = 0; ii < shift/frame_shift_samp;){
      if (AddParameters(&f0[ring_buffer_index+ii+cut], 1,
        &spectrogram[ring_buffer_index+ii+cut], &aperiodicity[ring_buffer_index+ii+cut],
        &synthesizer) == 1) {++ii;}

      while (Synthesis2(&synthesizer) != 0) {
        for (int j = 0; j < buffer_size; ++j){
          y[j + index] = synthesizer.buffer[j];
        }
        index += buffer_size;
      }

      if (IsLocked(&synthesizer) == 1) {
        printf("Locked!\n");
        break;
      }
    }

    ring_buffer_index += f0_length;
    if(ring_buffer_index>=ring_buffer_size) ring_buffer_index = 0;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  cout << "elapsed: " << elapsed/1000.0 << "[sec]" << endl;
}

void doit(double *x, int x_length, int fs, double *y){
  HarvestOption h_option = { 0 };
  InitializeHarvestOption(&h_option);
  h_option.frame_period = 4.0;
  h_option.f0_floor = 80.0;

  //Dio
  DioOption dio_option = {0};
  InitializeDioOption(&dio_option);
  dio_option.frame_period = 4.0;
  dio_option.speed = 1;
  dio_option.f0_floor = 40.0;
  dio_option.allowed_range = 0.1;

  int overlap = 2048;
  int shift=2048;
  int frame_shift_samp = 4*fs/1000;
  int win_length=shift+overlap;
  // int f0_length = GetSamplesForHarvest(fs, x_length, h_option.frame_period);
  int f0_length = GetSamplesForDIO(fs, x_length, h_option.frame_period);
  cout << f0_length << ", " << x_length << endl;
  // int real_f0_length = GetSamplesForHarvest(fs, win_length, h_option.frame_period);
  int real_f0_length = GetSamplesForDIO(fs, win_length, h_option.frame_period);
  cout << "real: " << real_f0_length << endl;

  double *f0 = new double[f0_length];
  double *refined_f0 = new double[f0_length];
  double *time_axis = new double[f0_length];

  double *real_f0=new double[f0_length+real_f0_length];
  double *tmp_f0 = new double[2*real_f0_length];
  double *pre_f0 = new double[2*real_f0_length];
  CheapTrickOption c_option = {0};
  InitializeCheapTrickOption(fs, &c_option);
  c_option.f0_floor = 71.0;
  c_option.fft_size = GetFFTSizeForCheapTrick(fs, &c_option);
  cout << "Cheap Trick FFT length: "<< c_option.fft_size<<endl;
  double **spectrogram = new double *[f0_length];
  double **real_spectrogram = new double *[f0_length+real_f0_length];
  double **tmp_spectrogram = new double *[2*real_f0_length];
  int hfft = c_option.fft_size/2+1;
  for(int i=0; i<f0_length; i++){
    spectrogram[i] = new double[hfft];
  }
  for(int i=0; i<f0_length+real_f0_length; i++){
    real_spectrogram[i] = new double[hfft];
  }
  for(int i=0; i<2*real_f0_length; i++){
    tmp_spectrogram[i] = new double[hfft];
  }
  D4COption d_option = {0};
  InitializeD4COption(&d_option);
  d_option.threshold = 0.85;
  double **aperiodicity = new double *[f0_length];
  double **real_aperiodicity = new double *[f0_length+real_f0_length];
  double **tmp_aperiodicity = new double *[2*real_f0_length];
  for(int i=0; i<f0_length; i++){
    aperiodicity[i] = new double[hfft];
  }
  for(int i=0; i<f0_length+real_f0_length; i++){
    real_aperiodicity[i] = new double[hfft];
  }
  for(int i=0; i<2*real_f0_length; i++){
    tmp_aperiodicity[i] = new double[hfft];
  }

  WorldSynthesizer synthesizer = { 0 };
  int buffer_size = 256;
  InitializeSynthesizer(fs, h_option.frame_period,
  c_option.fft_size, buffer_size, 200, &synthesizer);

  printf("\nAnalysis\n");
  // Harvest(x, x_length, fs, &h_option, time_axis, f0);
  Dio(x, x_length, fs, &dio_option, time_axis, f0);
  StoneMask(x, x_length, fs, time_axis,
      f0, f0_length, refined_f0);
  for(int i=0;i<f0_length;i++){
    f0[i]=refined_f0[i];
  }
  CheapTrick(x, x_length, fs, time_axis, f0, f0_length, &c_option, spectrogram);
  D4C(x, x_length, fs, time_axis, f0, f0_length, c_option.fft_size, &d_option, aperiodicity);
  cout<<endl;
  printf("\nAnalysis\n");
  int cut = overlap/2/frame_shift_samp;
  cout << cut << endl;
  cout << frame_shift_samp << endl;
  cout << real_f0_length << endl;

  int index = 0;
  for(int i=0; i<x_length-win_length;i+=shift){
    // F0
    // Harvest(&x[i], win_length, fs, &h_option, time_axis, tmp_f0);
    Dio(&x[i], win_length, fs, &dio_option, time_axis, pre_f0);
    // cout << "dio" << endl;
    StoneMask(&x[i], win_length, fs, time_axis,
      pre_f0, f0_length, tmp_f0);
    // cout << "dio2" << endl;
    for(int j=0; j<shift/frame_shift_samp; j++){
      real_f0[i/frame_shift_samp+j+cut]=tmp_f0[j+cut];
      // cout << tmp_f0[j+cut] << ",";
    }
    // cout << endl;
    if(i==0){
      for(int j=0; j<cut; j++){
        real_f0[j]=tmp_f0[j];
      }
    }
    // cout << "sp" << endl;

    // Sp
    // SpectralEnvelopeEstimation(x, x_length, &world_parameters);
    CheapTrick(&x[i], win_length, fs, time_axis, tmp_f0, real_f0_length, &c_option, tmp_spectrogram);
    // cout << "sp" << endl;
    for(int j=0; j<shift/frame_shift_samp; j++){
      // cout << j << ",";
      for(int k=0;k<hfft; k++){
        real_spectrogram[i/frame_shift_samp+j+cut][k]=tmp_spectrogram[j+cut][k];
      }
    }
    if(i==0){
      for(int j=0; j<cut; j++){
        for(int k=0;k<hfft; k++){
          real_spectrogram[j][k]=tmp_spectrogram[j][k];
        }
      }
    }
    
    // Ap
    // AperiodicityEstimation(x, x_length, &world_parameters);
    D4C(x, win_length, fs, time_axis, tmp_f0, real_f0_length, c_option.fft_size, &d_option, tmp_aperiodicity);
    for(int j=0; j<shift/frame_shift_samp; j++){
      for(int k=0;k<hfft; k++){
        real_aperiodicity[i/frame_shift_samp+j+cut][k]=tmp_aperiodicity[j+cut][k];
      }
    }
    if(i==0){
      for(int j=0; j<cut; j++){
        for(int k=0;k<hfft; k++){
          real_aperiodicity[j][k]=tmp_aperiodicity[j][k];
        }
      }
    }

    // synth
    
    for (int j = 0; j < shift/frame_shift_samp;) {
      // works
      // Add one frame (i shows the frame index that should be added)
      if (AddParameters(&real_f0[i/frame_shift_samp+j+cut], 1,
        &real_spectrogram[i/frame_shift_samp+j+cut], &real_aperiodicity[i/frame_shift_samp+j+cut],
        &synthesizer) == 1) {++j;}

      // not works
      // tmp_f0[j+cut]=real_f0[i/frame_shift_samp+j+cut];
      // for(int k=0;k<hfft;k++){
      //   tmp_spectrogram[j+cut][k]=real_spectrogram[i/frame_shift_samp+j+cut][k];
      //   tmp_aperiodicity[j+cut][k]=real_aperiodicity[i/frame_shift_samp+j+cut][k];
      // }
      // if (AddParameters(&tmp_f0[j+cut], 1,
      //   &tmp_spectrogram[j+cut], &tmp_aperiodicity[j+cut],
      //   &synthesizer) == 1) {++j;}

      // Synthesize speech with length of buffer_size sample.
      // It is repeated until the function returns 0
      // (it suggests that the synthesizer cannot generate speech).
      while (Synthesis2(&synthesizer) != 0) {
        for (int j = 0; j < buffer_size; ++j){
          y[j + index] = synthesizer.buffer[j];
        }
        index += buffer_size;
      }

      // Check the "Lock" (Please see synthesisrealtime.h)
      if (IsLocked(&synthesizer) == 1) {
        printf("Locked!\n");
        break;
      }
    }
  }
  // for(int i=0; i< f0_length; i++){
  //   cout<<real_f0[i]<<",";
  // }
  cout<<endl;
  // CheapTrick(x, x_length, fs, time_axis, real_f0, f0_length, &c_option, real_spectrogram);

  cout<<endl;
  // for (int ii = 0; ii < f0_length;) {
  //   // Add one frame (i shows the frame index that should be added)
  //   if (AddParameters(&real_f0[ii], 1,
  //     &real_spectrogram[ii], &real_aperiodicity[ii],
  //     &synthesizer) == 1) {++ii;}

  //   // Synthesize speech with length of buffer_size sample.
  //   // It is repeated until the function returns 0
  //   // (it suggests that the synthesizer cannot generate speech).
  //   while (Synthesis2(&synthesizer) != 0) {
  //     for (int j = 0; j < buffer_size; ++j){
  //       y[j + index] = synthesizer.buffer[j];
  //     }
  //     index += buffer_size;
  //   }

  //   // Check the "Lock" (Please see synthesisrealtime.h)
  //   if (IsLocked(&synthesizer) == 1) {
  //     printf("Locked!\n");
  //     break;
  //   }
  // }

  // f0 compair
  for(int i=0; i< f0_length; i++){
    // cout << (f0[i]-real_f0[i])*(f0[i]-real_f0[i])<<","<<endl;
    cout<<"("<<f0[i]<<","<<real_f0[i]<<"),";
  }
  cout<<endl;
  for(int i=0; i< f0_length; i++){
    //sq err
    double sum = 0.0;
    double psum = 0.0;
    // cout << "j: " << i << ", ";
    for (int k=0; k<hfft; k++){
      // cout << "k: " << k << ", ";
      // sum+=(real_spectrogram[i][k]-spectrogram[i][k])*(real_spectrogram[i][k]-spectrogram[i][k])/(spectrogram[i][k]*spectrogram[i][k]+1e-20);
      // sum+=(real_spectrogram[i][k]-spectrogram[i][k])*(real_spectrogram[i][k]-spectrogram[i][k]);
      sum+=real_spectrogram[i][k]*real_spectrogram[i][k];
      psum+=spectrogram[i][k]*spectrogram[i][k];
      // sum+=real_aperiodicity[i][k]*real_aperiodicity[i][k];
      // psum+=aperiodicity[i][k]*aperiodicity[i][k];
    }
    cout << "(" << sum << ", " << psum << ")";
    // cout <<  ", " << psum ;
  }
  cout<<endl;
}


namespace {

void DisplayInformation(int fs, int nbit, int x_length) {
  printf("File information\n");
  printf("Sampling : %d Hz %d Bit\n", fs, nbit);
  printf("Length %d [sample]\n", x_length);
  printf("Length %f [sec]\n", static_cast<double>(x_length) / fs);
}

void F0EstimationDio(double *x, int x_length,
    WorldParameters *world_parameters) {
  DioOption option = {0};
  InitializeDioOption(&option);

  // Modification of the option
  option.frame_period = world_parameters->frame_period;

  // Valuable option.speed represents the ratio for downsampling.
  // The signal is downsampled to fs / speed Hz.
  // If you want to obtain the accurate result, speed should be set to 1.
  option.speed = 1;

  // You can set the f0_floor below world::kFloorF0.
  option.f0_floor = 40.0;

  // You can give a positive real number as the threshold.
  // Most strict value is 0, but almost all results are counted as unvoiced.
  // The value from 0.02 to 0.2 would be reasonable.
  option.allowed_range = 0.1;

  // Parameters setting and memory allocation.
  world_parameters->f0_length = GetSamplesForDIO(world_parameters->fs,
    x_length, world_parameters->frame_period);
  world_parameters->f0 = new double[world_parameters->f0_length];
  world_parameters->time_axis = new double[world_parameters->f0_length];
  double *refined_f0 = new double[world_parameters->f0_length];

  printf("\nAnalysis\n");
  Dio(x, x_length, world_parameters->fs, &option, world_parameters->time_axis,
      world_parameters->f0);

  // StoneMask is carried out to improve the estimation performance.
  StoneMask(x, x_length, world_parameters->fs, world_parameters->time_axis,
      world_parameters->f0, world_parameters->f0_length, refined_f0);

  for (int i = 0; i < world_parameters->f0_length; ++i)
    world_parameters->f0[i] = refined_f0[i];

  delete[] refined_f0;
}

// void F0EstimationHarvest(double *x, int x_length,
//     WorldParameters *world_parameters) {
//   HarvestOption option = { 0 };
//   InitializeHarvestOption(&option);

//   // You can change the frame period.
//   // But the estimation is carried out with 1-ms frame shift.
//   option.frame_period = world_parameters->frame_period;

//   // You can set the f0_floor below world::kFloorF0.
//   option.f0_floor = 40.0;

//   // Parameters setting and memory allocation.
//   world_parameters->f0_length = GetSamplesForHarvest(world_parameters->fs,
//     x_length, world_parameters->frame_period);
//   world_parameters->f0 = new double[world_parameters->f0_length];
//   world_parameters->time_axis = new double[world_parameters->f0_length];

//   printf("\nAnalysis\n");
//   Harvest(x, x_length, world_parameters->fs, &option,
//       world_parameters->time_axis, world_parameters->f0);
// }

// void SpectralEnvelopeEstimation(double *x, int x_length,
//     WorldParameters *world_parameters) {
//   CheapTrickOption option = {0};
//   // Note (2017/01/02): fs is added as an argument.
//   InitializeCheapTrickOption(world_parameters->fs, &option);

//   // Default value was modified to -0.15.
//   // option.q1 = -0.15;

//   // Important notice (2017/01/02)
//   // You can set the fft_size.
//   // Default is GetFFTSizeForCheapTrick(world_parameters->fs, &option);
//   // When fft_size changes from default value,
//   // a replaced f0_floor will be used in CheapTrick().
//   // The lowest F0 that WORLD can work as expected is determined
//   // by the following : 3.0 * fs / fft_size.
//   option.f0_floor = 71.0;
//   option.fft_size = GetFFTSizeForCheapTrick(world_parameters->fs, &option);
//   // We can directly set fft_size.
// //   option.fft_size = 1024;

//   // Parameters setting and memory allocation.
//   world_parameters->fft_size = option.fft_size;
//   world_parameters->spectrogram = new double *[world_parameters->f0_length];
//   for (int i = 0; i < world_parameters->f0_length; ++i)
//     world_parameters->spectrogram[i] =
//       new double[world_parameters->fft_size / 2 + 1];

//   CheapTrick(x, x_length, world_parameters->fs, world_parameters->time_axis,
//       world_parameters->f0, world_parameters->f0_length, &option,
//       world_parameters->spectrogram);
// }

// void AperiodicityEstimation(double *x, int x_length,
//     WorldParameters *world_parameters) {
//   D4COption option = {0};
//   InitializeD4COption(&option);

//   // Comment was modified because it was confusing (2017/12/10).
//   // It is used to determine the aperiodicity in whole frequency band.
//   // D4C identifies whether the frame is voiced segment even if it had an F0.
//   // If the estimated value falls below the threshold,
//   // the aperiodicity in whole frequency band will set to 1.0.
//   // If you want to use the conventional D4C, please set the threshold to 0.0.
//   option.threshold = 0.85;

//   // Parameters setting and memory allocation.
//   world_parameters->aperiodicity = new double *[world_parameters->f0_length];
//   for (int i = 0; i < world_parameters->f0_length; ++i)
//     world_parameters->aperiodicity[i] =
//       new double[world_parameters->fft_size / 2 + 1];

//   D4C(x, x_length, world_parameters->fs, world_parameters->time_axis,
//       world_parameters->f0, world_parameters->f0_length,
//       world_parameters->fft_size, &option, world_parameters->aperiodicity);
// }

// void ParameterModification(int argc, char *argv[], int fs, int f0_length,
//     int fft_size, double *f0, double **spectrogram) {
//   // F0 scaling
//   if (argc >= 4) {
//     double shift = atof(argv[3]);
//     for (int i = 0; i < f0_length; ++i) f0[i] *= shift;
//   }
//   if (argc < 5) return;

//   // Spectral stretching
//   double ratio = atof(argv[4]);
//   double *freq_axis1 = new double[fft_size];
//   double *freq_axis2 = new double[fft_size];
//   double *spectrum1 = new double[fft_size];
//   double *spectrum2 = new double[fft_size];

//   for (int i = 0; i <= fft_size / 2; ++i) {
//     freq_axis1[i] = ratio * i / fft_size * fs;
//     freq_axis2[i] = static_cast<double>(i) / fft_size * fs;
//   }
//   for (int i = 0; i < f0_length; ++i) {
//     for (int j = 0; j <= fft_size / 2; ++j)
//       spectrum1[j] = log(spectrogram[i][j]);
//     interp1(freq_axis1, spectrum1, fft_size / 2 + 1, freq_axis2,
//       fft_size / 2 + 1, spectrum2);
//     for (int j = 0; j <= fft_size / 2; ++j)
//       spectrogram[i][j] = exp(spectrum2[j]);
//     if (ratio >= 1.0) continue;
//     for (int j = static_cast<int>(fft_size / 2.0 * ratio);
//         j <= fft_size / 2; ++j)
//       spectrogram[i][j] =
//       spectrogram[i][static_cast<int>(fft_size / 2.0 * ratio) - 1];
//   }
//   delete[] spectrum1;
//   delete[] spectrum2;
//   delete[] freq_axis1;
//   delete[] freq_axis2;
// }

// void WaveformSynthesis(WorldParameters *world_parameters, int fs,
//     int y_length, double *y) {
//   // Synthesis by the aperiodicity
//   printf("\nSynthesis 1 (conventional algorithm)\n");
//   Synthesis(world_parameters->f0, world_parameters->f0_length,
//       world_parameters->spectrogram, world_parameters->aperiodicity,
//       world_parameters->fft_size, world_parameters->frame_period, fs,
//       y_length, y);
// }

// void WaveformSynthesis2(WorldParameters *world_parameters, int fs,
//     int y_length, double *y) {
//   printf("\nSynthesis 2 (All frames are added at the same time)\n");

//   WorldSynthesizer synthesizer = { 0 };
//   int buffer_size = 64;
//   InitializeSynthesizer(world_parameters->fs, world_parameters->frame_period,
//       world_parameters->fft_size, buffer_size, 1, &synthesizer);

//   // All parameters are added at the same time.
//   AddParameters(world_parameters->f0, world_parameters->f0_length,
//       world_parameters->spectrogram, world_parameters->aperiodicity,
//       &synthesizer);

//   int index;
//   for (int i = 0; Synthesis2(&synthesizer) != 0; ++i) {
//     index = i * buffer_size;
//     for (int j = 0; j < buffer_size; ++j)
//       y[j + index] = synthesizer.buffer[j];
//   }

//   DestroySynthesizer(&synthesizer);
// }

// void WaveformSynthesis3(WorldParameters *world_parameters, int fs,
//     int y_length, double *y) {
//   // Synthesis by the aperiodicity
//   printf("\nSynthesis 3 (Ring buffer is efficiently used.)\n");

//   WorldSynthesizer synthesizer = { 0 };
//   int buffer_size = 64;
//   InitializeSynthesizer(world_parameters->fs, world_parameters->frame_period,
//       world_parameters->fft_size, buffer_size, 100, &synthesizer);

//   int offset = 0;
//   int index = 0;
//   for (int i = 0; i < world_parameters->f0_length;) {
//     // Add one frame (i shows the frame index that should be added)
//     if (AddParameters(&world_parameters->f0[i], 1,
//       &world_parameters->spectrogram[i], &world_parameters->aperiodicity[i],
//       &synthesizer) == 1) ++i;

//     // Synthesize speech with length of buffer_size sample.
//     // It is repeated until the function returns 0
//     // (it suggests that the synthesizer cannot generate speech).
//     while (Synthesis2(&synthesizer) != 0) {
//       index = offset * buffer_size;
//       for (int j = 0; j < buffer_size; ++j)
//         y[j + index] = synthesizer.buffer[j];
//       offset++;
//     }

//     // Check the "Lock" (Please see synthesisrealtime.h)
//     if (IsLocked(&synthesizer) == 1) {
//       printf("Locked!\n");
//       break;
//     }
//   }

//   DestroySynthesizer(&synthesizer);
// }

// void DestroyMemory(WorldParameters *world_parameters) {
//   delete[] world_parameters->time_axis;
//   delete[] world_parameters->f0;
//   for (int i = 0; i < world_parameters->f0_length; ++i) {
//     delete[] world_parameters->spectrogram[i];
//     delete[] world_parameters->aperiodicity[i];
//   }
//   delete[] world_parameters->spectrogram;
//   delete[] world_parameters->aperiodicity;
// }

}  // namespace

//-----------------------------------------------------------------------------
// Test program.
// test.exe input.wav outout.wav f0 spec flag
// input.wav  : argv[1] Input file
// output.wav : argv[2] Output file
// f0         : argv[3] F0 scaling (a positive number)
// spec       : argv[4] Formant shift (a positive number)
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3 && argc != 4 && argc != 5) {
    printf("error\n");
    return -2;
  }

  // Memory allocation is carried out in advanse.
  // This is for compatibility with C language.
  int x_length = GetAudioLength(argv[1]);
  if (x_length <= 0) {
    if (x_length == 0) printf("error: File not found.\n");
    else printf("error: The file is not .wav format.\n");
    return -1;
  }
  double *x = new double[x_length];
  // wavread() must be called after GetAudioLength().
  int fs, nbit;
  wavread(argv[1], &fs, &nbit, x);
  DisplayInformation(fs, nbit, x_length);

  // //---------------------------------------------------------------------------
  // // Analysis part
  // //---------------------------------------------------------------------------
  // WorldParameters world_parameters = { 0 };
  // // You must set fs and frame_period before analysis/synthesis.
  // world_parameters.fs = fs;
  // // 5.0 ms is the default value.
  // world_parameters.frame_period = 5.0;

  // // F0 estimation
  // // DIO
  // F0EstimationDio(x, x_length, &world_parameters);

  // // Harvest
  // // F0EstimationHarvest(x, x_length, &world_parameters);

  // // Spectral envelope estimation
  // SpectralEnvelopeEstimation(x, x_length, &world_parameters);

  // // Aperiodicity estimation by D4C
  // AperiodicityEstimation(x, x_length, &world_parameters);

  // // Note that F0 must not be changed until all parameters are estimated.
  // ParameterModification(argc, argv, fs, world_parameters.f0_length,
  //     world_parameters.fft_size, world_parameters.f0,
  //     world_parameters.spectrogram);

  // //---------------------------------------------------------------------------
  // // Synthesis part
  // // There are three samples in speech synthesis
  // // 1: Conventional synthesis
  // // 2: Example of real-time synthesis
  // // 3: Example of real-time synthesis (Ring buffer is efficiently used)
  // //---------------------------------------------------------------------------
  char filename[100];
  // // The length of the output waveform
  // int y_length = static_cast<int>((world_parameters.f0_length - 1) *
  //   world_parameters.frame_period / 1000.0 * fs) + 1;
  // double *y = new double[y_length];

  // Synthesis 1 (conventional synthesis)
  // for (int i = 0; i < y_length; ++i) y[i] = 0.0;
  // WaveformSynthesis(&world_parameters, fs, y_length, y);
  // sprintf(filename, "01%s", argv[2]);
  // wavwrite(y, y_length, fs, 16, filename);

  // // Synthesis 2 (All frames are added at the same time)
  // for (int i = 0; i < y_length; ++i) y[i] = 0.0;
  // WaveformSynthesis2(&world_parameters, fs, y_length, y);
  // sprintf(filename, "02%s", argv[2]);
  // wavwrite(y, y_length, fs, 16, filename);

  // Synthesis 3 (Ring buffer is efficiently used.)
  // for (int i = 0; i < y_length; ++i) y[i] = 0.0;
  // WaveformSynthesis3(&world_parameters, fs, y_length, y);
  // sprintf(filename, "%s", argv[2]);
  // wavwrite(y, y_length, fs, 16, filename);

  // mine
  double* y = new double[x_length];
  for(int k=0;k<x_length;k++){
    y[k]=0.0;
  }
  // doit(x, x_length, fs, y);
  online(x, x_length, fs, y);
  sprintf(filename, "%s", argv[2]);
  wavwrite(y, x_length, fs, 16, filename);

  // delete[] y;
  // delete[] x;
  // DestroyMemory(&world_parameters);

  printf("complete.\n");
  return 0;
}
