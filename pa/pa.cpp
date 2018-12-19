#include <iostream>
#include <chrono>

#include "portaudio.h"

#include "world/d4c.h"
#include "world/harvest.h"
// #include "world/matlabfunctions.h"
#include "world/cheaptrick.h"
#include "world/synthesisrealtime.h"

#include<stdio.h>
#include<math.h>
#include"portaudio.h"
#define Fs 16000 //サンプリング周波数
#define FRAMES_PER_BUFFER 2048 //バッファサイズ

/*ユーザ定義データ*/
typedef struct{
    float freq; //正弦波の周波数
    float index;
}padata;

struct WorldOptions{
  int fs;
  int overlap;
  int shift;
  int frame_shift_samp;
  int win_length;
  int f0_length;
  int hfft;
  int cut;
  int buffer_size;
  int ring_buffer_size;
  int fft_size;
  double frame_period;
  double f0_floor;
  double c_f0_floor;
  HarvestOption h_option;
  CheapTrickOption c_option;
  D4COption d_option;

  WorldOptions(){}
  WorldOptions(int fs);
  void defaultOptions();
};

struct WorldParams{
  int ring_buffer_index;
  WorldOptions options;
  
  double *buffer;
  double *f0;
  double *time_axis;
  double **spectrogram;
  double **aperiodicity;
  WorldSynthesizer synthesizer;
  WorldParams();
  void initParams(const WorldOptions &options);
};

using namespace std;

WorldOptions::WorldOptions(int fs){
  this->fs = fs;
  this->overlap = 2048;
  this->shift=2048;
  this->win_length=shift+overlap;
  this->frame_period = 4.0;
  this->frame_shift_samp = static_cast<int>(frame_period*fs/1000);
  this->cut = overlap/2/frame_shift_samp;
  this->f0_length = GetSamplesForHarvest(fs, win_length, frame_period);
  this->buffer_size = 64;
  this->ring_buffer_size = 4*f0_length;
  this->f0_floor = 80.0;
  this->c_f0_floor = 71.0;
  this->h_option = {0};
  this->c_option = {0};
  this->d_option = {0};
}

void WorldOptions::defaultOptions(){
  InitializeHarvestOption(&this->h_option);
  this->h_option.frame_period = this->frame_period;
  this->h_option.f0_floor = this->f0_floor;
  InitializeCheapTrickOption(fs, &this->c_option);
  this->c_option.f0_floor = this->c_f0_floor;
  this->fft_size = GetFFTSizeForCheapTrick(fs, &this->c_option);
  hfft = fft_size/2+1;
  this->c_option.fft_size = this->fft_size;
  InitializeD4COption(&this->d_option);
  this->d_option.threshold = 0.85;
}

WorldParams::WorldParams(){
  synthesizer = {0};
}

void WorldParams::initParams(const WorldOptions &options){
  ring_buffer_index = 0;
  this->options = options;

  buffer = new double[options.win_length];
  for(int i=0;i<options.win_length;i++){
    buffer[i]=0.0;
  }
  f0 = new double[options.ring_buffer_size];
  time_axis = new double[options.f0_length];
  spectrogram = new double *[options.ring_buffer_size];
  aperiodicity = new double *[options.ring_buffer_size];
  for(int i=0; i<options.ring_buffer_size; i++){
    f0[i]=0.0;
    spectrogram[i] = new double[options.hfft];
    aperiodicity[i] = new double[options.hfft];
    for(int j=0;j<options.hfft;j++){
      spectrogram[i][j]=1e-16;
      aperiodicity[i][j]=1e-16;
    }
  }
  InitializeSynthesizer(options.fs, options.frame_period, options.fft_size, options.buffer_size, options.ring_buffer_size, &synthesizer);
}

void online(const double *x, int x_length, int fs, double *y){
  double frame_period = 4.0;
  double f0_floor = 80.0;
  double c_f0_floor = 71.0;
  HarvestOption h_option = { 0 };
  CheapTrickOption c_option = {0};
  D4COption d_option = {0};
  InitializeHarvestOption(&h_option);
  h_option.frame_period = frame_period;
  h_option.f0_floor = f0_floor;
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
  double *time_axis = new double[f0_length];
  double **spectrogram = new double *[ring_buffer_size];
  double **aperiodicity = new double *[ring_buffer_size];
  for(int i=0; i<ring_buffer_size; i++){
    spectrogram[i] = new double[hfft];
    aperiodicity[i] = new double[hfft];
  }
  
  InitializeSynthesizer(fs, frame_period, fft_size, buffer_size, ring_buffer_size, &synthesizer);

  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
  int index = 0;
  int ring_buffer_index = 0;
  for(int i=0; i<x_length-win_length;i+=shift){
    Harvest(&x[i], win_length, fs, &h_option, time_axis, &f0[ring_buffer_index]);
    CheapTrick(&x[i], win_length, fs, time_axis, &f0[ring_buffer_index], f0_length, &c_option, &spectrogram[ring_buffer_index]);
    D4C(&x[i], win_length, fs, time_axis, &f0[ring_buffer_index], f0_length, fft_size, &d_option, &aperiodicity[ring_buffer_index]);

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

void convert(double f0mul, double spmul, double* f0, double** sp, double** ap, int f0_length, int hfft){
  for(int i=0; i<f0_length; i++){
    f0[i]*=f0mul;
    cout << f0[i] << endl;
  }
}


/* オーディオ処理コールバック関数*/
static int dsp(const void *inputBuffer, //入力
               void *outputBuffer, //出力
               unsigned long framesPerBuffer,
               const PaStreamCallbackTimeInfo *timeInfo,
               PaStreamCallbackFlags statusFlags,
               void *userData //ユーザ定義データ 
               ){
  // std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
  WorldParams *data = (WorldParams *)userData;
  float *in = (float *)inputBuffer;
  float *out = (float *)outputBuffer;
  long i, j, ii, index;
  int shift = data->options.shift;
  int overlap = data->options.overlap;
  int win_length = data->options.win_length;
  double* buffer = data->buffer;
  int fs = data->options.fs;
  double *time_axis = data->time_axis;
  double *f0 = data->f0;
  int ring_buffer_index = data->ring_buffer_index;
  int f0_length = data->options.f0_length;
  double** spectrogram = data->spectrogram;
  double** aperiodicity = data->aperiodicity;
  int fft_size = data->options.fft_size;
  int frame_shift_samp = data->options.frame_shift_samp;
  int cut = data->options.cut;
  int buffer_size = data->options.buffer_size;

  // for( i=0; i<framesPerBuffer; i++){
  //     out[i] = in[i];
  // }
  // return 0;

  // shift and copy
  for(i=0; i<overlap; i++){
    buffer[i] = buffer[i+shift];
  }
  for(i=0;i<shift;i++){
    buffer[i+shift] = in[i];
  }

  // analysis
  Harvest(buffer, win_length, fs, &data->options.h_option, time_axis, &f0[ring_buffer_index]);
  
  CheapTrick(buffer, win_length, fs, time_axis, &f0[ring_buffer_index], f0_length, &data->options.c_option, &spectrogram[ring_buffer_index]);
  D4C(buffer, win_length, fs, time_axis, &f0[ring_buffer_index], f0_length, fft_size, &data->options.d_option, &aperiodicity[ring_buffer_index]);

  // convert
  convert(2.0, 1.2, &f0[ring_buffer_index], spectrogram, aperiodicity, f0_length, data->options.hfft);

  // synthesis
  index = 0;
  for (ii = 0; ii < shift/frame_shift_samp;){
    if (AddParameters(&f0[ring_buffer_index+ii+cut], 1,
      &spectrogram[ring_buffer_index+ii+cut], &aperiodicity[ring_buffer_index+ii+cut],
      &data->synthesizer) == 1) {++ii;}

    while (Synthesis2(&data->synthesizer) != 0) {
      for (j = 0; j < buffer_size; ++j){
        out[j + index] = data->synthesizer.buffer[j];
      }
      index += buffer_size;
    }

    if (IsLocked(&data->synthesizer) == 1) {
      printf("Locked!\n");
      break;
    }
  }

  data->ring_buffer_index += data->options.f0_length;
  if(data->ring_buffer_index >= data->options.ring_buffer_size) data->ring_buffer_index = 0;

  // std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  // double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  // cout << "elapsed: " << elapsed/1000.0 << "[sec]" << endl;
  return 0;
}

int main(){
  constexpr int fs = 16000;
    PaStreamParameters inParam;
    PaStreamParameters outParam;
    PaStream *stream;
    PaError err;
    WorldOptions options(fs);
    options.defaultOptions();
    WorldParams data;
    data.initParams(options);

    //PortAudio初期化
    Pa_Initialize();

    inParam.device = Pa_GetDefaultInputDevice();
    inParam.channelCount = 1;
    inParam.sampleFormat = paFloat32;
    inParam.suggestedLatency = Pa_GetDeviceInfo( inParam.device )->defaultLowOutputLatency;
    inParam.hostApiSpecificStreamInfo = nullptr;
    //出力の設定
    outParam.device = Pa_GetDefaultOutputDevice();
    outParam.channelCount = 1;
    outParam.sampleFormat = paFloat32;
    outParam.suggestedLatency = Pa_GetDeviceInfo( outParam.device )->defaultLowOutputLatency;
    outParam.hostApiSpecificStreamInfo = nullptr;

    //PortAudioオープン
    Pa_OpenStream(
        &stream,
        &inParam,
        &outParam,
        Fs,
        FRAMES_PER_BUFFER,
        paClipOff,
        dsp,
        &data);
    
    //PortAudioスタート
    Pa_StartStream(stream);

    //エンターキーが押されるまで待機
    getchar();

    //PortAudio終了
    Pa_StopStream(stream);
    Pa_CloseStream(stream);
    Pa_Terminate();
    return 0;
    return 0;
}