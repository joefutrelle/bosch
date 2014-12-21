#include <stdio.h>

#include <stdlib.h>
#include <math.h>

#include <tiffio.h>
#include <sndfile.h>

#include <pthread.h>

#define TABLE_SIZE 1024
#define PI 3.1415926535
#define SAMPLE_RATE 44100

#define MODE_COAGULA 1
#define MODE_BISCH 1

#define VERSION "1.2"
#define VERSION_DATE "May 2009"

int coagulaMode = 0;
int bischMode = 0;

/* transition time between volume levels */
/* #define RAMP_TIME 441 */
/* 5ms (too fast?) */
#define RAMP_TIME 220

#define SAMPLE float
#define SFWRITE sf_write_float

#define N_CHANNELS 2

#define MAX_THREADS 16

/* number of different waveforms */
#define N_WAVE_TYPES 3

#define WAVE_SINE 0
#define WAVE_SQUARE 1
#define WAVE_SAW 2

#define WAVE_RAND 2

int nThreads = 8; /* global: number of threads to use */

float wav_tab[N_WAVE_TYPES][TABLE_SIZE];
/* below for temporary */
float ***wav_tab_scaled;

/**
 * Converts hz to wavelength given the
 * sampling rate
 */
#define HZ2WLFORSL(hz,sr) ((int)((float)sr / (float)hz))
#define HZ2WL(hz) ((int)(44100.0 / (float)hz))

#define BLOCK_SIZE 32768
SAMPLE buffer[BLOCK_SIZE * N_CHANNELS];
int buf_ix = 0;

#define ENERGY_THRESHOLD 0.00005

void make_tables(int wave_type) {
     int i, wt;
     float x, step = (float)(2 * PI) / (float) TABLE_SIZE;
     int half = TABLE_SIZE / 2;
     float square_amp = 0.5;
     float rand_amp = 0.75;
     float prev_rand = 0.0;
     for(i = 0, x = 0; i < TABLE_SIZE; i++, x += step) {
	  wav_tab[WAVE_SINE]  [i] = sin(x);
	  if(!coagulaMode) {
	       wav_tab[WAVE_SQUARE][i] = (i < half) ? (0 - square_amp) : square_amp;
/*
  wav_tab[WAVE_SAW]   [i] =
  (i < half) ?
  (((float)i / (float)half) * 2.0) - 1.0 :
  (((float)(TABLE_SIZE - i) / (float)half) * 2.0) - 1.0;
*/
	       wav_tab[WAVE_RAND]  [i] =
		    ((((float) ((rand() % 2000) - 1000)) / 2000.0) * rand_amp
		     + prev_rand) / 2.0;
	       prev_rand = wav_tab[WAVE_RAND][i];
	  } else {
	       wav_tab[WAVE_SQUARE][i] = wav_tab[WAVE_SINE][i];
	       wav_tab[WAVE_RAND][i] = wav_tab[WAVE_SINE][i];
	  }
     }
}

typedef struct {
  int fi; /* starting frequency band */
  int fn; /* ending frequency band */
  int *pitch; /* 2d array with wavelength for each band (in samples) */
  float **energy; /* 2d array with energy for each band (1.0 = max, 0.0 = min) and for each color r,g,b */
  float **energy_prev; /* 2d what the energy in the band was last frame */
  int *offset; /* per band, how many samples into the frame to wait before the transition to the new energy */
  float *pan; /* per-band panning */
  long ss; /* starting sample in absolute frames (for phase alignment) */
  int ns; /* number of samples to generate */
  SAMPLE *buf; /* output buffer */
} TIME_SLICE;

/**
 * f0 = starting frequency band
 * f1 = ending frequency band
 * pitch = array with wavelength for each band (in samples)
 * energy = 2d array with energy for each band (1.0 = max, 0.0 = min) and for
 * each color r,g,b
 * energy_prev = 2d what the energy in the band was last frame
 * offset = how many samples into the frame to wait before the transition to
 *  the new energy
 * pan = panning of each freq band
 * ss = starting sample number (to align phase)
 * ns = number of samples to generate
 */
void *do_timeSegment(void *vts) {
     TIME_SLICE *ts = (TIME_SLICE *)vts;
     long cs; /* current sample */
     long es = ts->ss + ts->ns; /* end sample */
     int i;
     int s; /* sample # within result segment */
     for(s = 0; s < ts->ns * N_CHANNELS; s++) {
	  ts->buf[s] = 0.0;
     }
     for(i = ts->fi; i < ts->fn; i++) { /* for each freq band */
	  int rgb;
	  int wl = ts->pitch[i];
	  float scale = (float)TABLE_SIZE / (float)wl;
	  int start_trans; /* beginning of power transition */
	  int end_trans; /* end of power transition */
	  start_trans = ts->offset[i] - RAMP_TIME;	
	  if(start_trans < 0) start_trans = 0;
	  end_trans = start_trans + RAMP_TIME;
	  if(end_trans >= ts->ns) end_trans = ts->ns;
	  for(rgb = 0; rgb < 3; rgb++) {
	       if(ts->energy_prev[rgb][i] > ENERGY_THRESHOLD ||
		  ts->energy[rgb][i] > ENERGY_THRESHOLD) {
		    float e = ts->energy_prev[rgb][i];
		    for(cs = ts->ss, s = 0; cs < es; cs++, s++) {
			 float amp;
			 if(s > start_trans && s < end_trans) {
			      e = ts->energy_prev[rgb][i] +
				   (((float)(s - start_trans) / (float)RAMP_TIME) *
				    (ts->energy[rgb][i] - ts->energy_prev[rgb][i]));
			 } else if(s >= end_trans) {
			      e = ts->energy[rgb][i];
			 }
			 /* if(s > ts->offset[i]) e = ts->energy[rgb][i]; */
			 /* compute amplitude */
			 amp = (wav_tab_scaled[rgb][i][cs % wl] * e);
			 /* now distribute amp across channels */
			 ts->buf[s*2]   += amp * ts->pan[i];
			 ts->buf[s*2+1] += amp * (1.0 - ts->pan[i]);
			 /* are we clipping? */
			 if((ts->buf[s*2] > 1.0 || ts->buf[s*2] < -1.0) ||
			    (ts->buf[s*2+1] > 1.0 || ts->buf[s*2+1] < -1.0)) {
			   fprintf(stderr," [CLIPPING]\n");
			   exit(-1);
			 }
		    }
	       }
	  }
     }
}

void output(SNDFILE *out, SAMPLE *buf, int n) {
     int i;
     for(i = 0; i < n * N_CHANNELS; i++) {
	  buffer[buf_ix++] = buf[i];
	  if(buf_ix==BLOCK_SIZE) {
	       SFWRITE (out, buffer, BLOCK_SIZE);
	       buf_ix = 0;
	  }
     }
}

void final_output(SNDFILE *out) {
     SFWRITE (out, buffer, buf_ix);
     sf_close(out);
}

void show_progress(int seg, int nseg) {
     int i, n;
     float percent = ((float)seg / (float)nseg) * 100;
     fprintf(stderr,"\rcompleted: [");
     n = percent / 3;
     for(i = 0; i < n; i++) {
	  fputc('*',stderr);
     }
     for(; i < 33; i++) {
	  fputc(' ',stderr);
     }
     fprintf(stderr,"] %.2f%%",percent);
     fflush(stderr);
}

void usage() {
     fprintf(stderr,"usage: bosch [options] tifffile\n");
     fprintf(stderr,"options:\n");
     fprintf(stderr,"   -d duration (s)\n");
     fprintf(stderr,"   -o output file\n");
     fprintf(stderr,"   -l lowest pitch (Hz)\n");
     fprintf(stderr,"   -r pitch range (semitones)\n");
     fprintf(stderr,"   -c coagula mode (sine waves only)\n");
     fprintf(stderr,"   -i bisch mode (energy = delta intensity)\n");
     fprintf(stderr,"   -a amp scale\n");
     fprintf(stderr,"   -f output format [wav|aiff] (default: wav)\n");
     fprintf(stderr,"   -t number of threads to use (default: 4)\n");
     fprintf(stderr,"example:\n");
     fprintf(stderr,"   bosch -d 30 -o foo.wav -l 27.5 -r 88 foo.tiff\n");
     fprintf(stderr,"default values:\n   duration = 60\n   lowest pitch = 27.5\n   pitch range = 110\n");
     exit(-1);
}

void version() {
     fprintf(stderr,"bosch %s (%s)\n",VERSION,VERSION_DATE);
     fprintf(stderr,"(not c) nonexistent software group\n");
     fprintf(stderr,"uses libtiff, libsndfile\n");
}
int main(int argc, char **argv) {
  int i, j, k, l;

     /** output buffering */
     SAMPLE *result;
     long cs, s;

     /* how long each column lasts */
     int block;

     /** control params */

     float minFreq = 27.5;
     float pitchRange = 110.0;

     float duration = 60.0; /* how long the image lasts in seconds */
     float rolloff = 2.25; /* log scaling of amplitude of freq bands */
     float minEnergy = 0.1; /* the amplitude scale of the highest pitch band */
     float bandwidth; /* width of each band (assume linear) */
     float band;
     float ampScale = 1.0;

     /** tiff file params */
     TIFF *tif;
     char *tif_fn = NULL;
     uint32 w, h;
     uint32 *image;
     int c;

     /* output file params */
     SNDFILE *out;
     char *out_fn = NULL;

     /* pitches of the freq bands */
     int *pitch;
     /* offset of amp transition for each band */
     int *offset;
     /* energy of the freq bands */
     float *energy[3];
     float *energy_prev[3];
     /* panning of each freq band */
     float *pan;
     /* amp scale of each band (higher bands are quieter) */
     float *freqResponse;
     float x, y; /* temp vars */
     int rgb, wt;
     /* file format */
     int fileFormat = SF_FORMAT_WAV | SF_FORMAT_PCM_16;

     /* timeslice data for inter-thread communication */
     TIME_SLICE ts[MAX_THREADS];

     /* threads */
     pthread_t threads[MAX_THREADS];
     pthread_attr_t thread_attr;

     /* first, read the command line arguments */
     for(i = 1; i < argc; i++) {
	  if(!strcmp(argv[i],"--help")) {
	       usage();
	  } else if(!strcmp(argv[i],"--version")) {
	       version();
	  } else if(!strcmp(argv[i],"-d")) {
	       /* duration in seconds */
	       if(!(++i < argc)) { usage(); }
	       duration = atof(argv[i]);
	       if(duration < 0) { usage(); }
	  } else if(!strcmp(argv[i],"-o")) {
	       /* output file */
	       if(!(++i < argc)) { usage(); }
	       out_fn = argv[i];
	  } else if(!strcmp(argv[i],"-t")) {
	       /* number of threads */
	       if(!(++i < argc)) { usage(); }
	       nThreads = (int) atoi(argv[i]);
	       if(nThreads < 1 || nThreads > MAX_THREADS) {
		 fprintf(stderr,"silly number of threads: max=%d\n",MAX_THREADS);
		 usage();
	       }
	  } else if(!strcmp(argv[i],"-l")) {
	       /* lowest pitch */
	       if(!(++i < argc)) { usage(); }
	       minFreq = (float) atof(argv[i]);
	  } else if(!strcmp(argv[i],"-r")) {
	       /* pitch range in semitones */
	       if(!(++i < argc)) { usage(); }
	       pitchRange = (float) atof(argv[i]);
	  } else if(!strcmp(argv[i],"-a")) {
	       /* amp scale */
	       if(!(++i < argc)) { usage(); }
	       ampScale = (float) atof(argv[i]);
	  } else if(!strcmp(argv[i],"-c")) {
	       coagulaMode = 1;
	  } else if(!strcmp(argv[i],"-i")) {
	       bischMode = 1;
	  } else if(!strcmp(argv[i],"-f")) {
	       if(!(++i < argc)) { usage(); }
	       if(!strcmp(argv[i],"wav")) {
		    fileFormat = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	       } else if(!strcmp(argv[i],"aiff")) {
		    fileFormat = SF_FORMAT_AIFF | SF_FORMAT_PCM_16;
	       }
	  } else {
	       if(tif_fn) { usage(); }
	       tif_fn = argv[i];
	  }
     }
     if(!tif_fn) {
	  usage();
     }

     fprintf(stderr, "output file: %s\n",out_fn ? out_fn : "(none)");
     fprintf(stderr, "duration (s): %f\n",duration);
     fprintf(stderr, "tiff file: %s\n",tif_fn);
     fprintf(stderr, "low pitch (Hz): %.3f\n",minFreq);
     fprintf(stderr, "pitch range (semitones): %.2f\n",pitchRange);

     /* first try to open the output file */
     if(out_fn) {
	  SF_INFO sfi;
	  sfi.samplerate = 44100;
	  sfi.channels = 2;
	  sfi.format = fileFormat;
	  out = sf_open(out_fn, SFM_WRITE, &sfi);
     } else {
	  SF_INFO sfi;
	  sfi.samplerate = 44100;
	  sfi.channels = 2;
	  sfi.format = fileFormat;
	  out = sf_open_fd(fileno(stdout), SFM_WRITE, &sfi, 1);
     }

     /* first arg is TIFF file. read TIFF file into raster */
     tif = TIFFOpen(tif_fn,"r");
     TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
     TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);

     image = (uint32*) _TIFFmalloc(w * h * sizeof (uint32));
     if (!image) {
	  fprintf(stderr,"error opening tiff file!");
	  exit(-1);
     }

     if(!TIFFReadRGBAImage(tif, w, h, image, 0)) {
	  fprintf(stderr,"error reading tiff image!");
	  exit(-1);
     }

     fprintf(stderr,"opened tiff image width=%d height=%d\n", w, h);

     /* make the wavetables */
     make_tables(WAVE_SAW);

     /* we know how wide the image is. map that to time */
     fprintf(stderr, "Each column is %.3f s\n",duration / (float)w);
     block = (duration / (float)w) * SAMPLE_RATE;
     result = calloc(block * N_CHANNELS * nThreads, sizeof(SAMPLE));

     /* ok, now we know how high the image is. map that to frequencies */
     pitch = calloc(h, sizeof(int));
     /* log scale so low bands are narrower */
     for(i = 0; i < h; i++) {
	  band = minFreq *
	       (pow(2.0, (((float)i/(float)h) * pitchRange) / 12.0));
	  pitch[i] = HZ2WL(band);
     }

     for(c = 0; c < 3; c++) {
	  energy[c] = calloc(h * N_WAVE_TYPES, sizeof(float));
	  energy_prev[c] = calloc(h * N_WAVE_TYPES, sizeof(float));
     }

     cs = 0;

     /* set random offsets to avoid "buzzing" */
     /* also set panning to random */
     /* also set up frequency response */
     offset = calloc(h, sizeof(int));
     pan = calloc(h, sizeof(float));
     freqResponse = calloc(h, sizeof(float));
     for(i = 0; i < h; i++) {
	  offset[i] = random() % block;
	  pan[i] = ((float) (random() % 1000)) / 1000.0;
	  freqResponse[i] =
	       (pow(1.0 - ((float)i / (float)h), rolloff) *
		(1.0 - minEnergy)) +
	       minEnergy;
     }

     /* now precompute all the wavetables */
     wav_tab_scaled = (float ***)calloc(N_WAVE_TYPES,sizeof(float **));
     for(rgb = 0; rgb < 3; rgb++) {
	  wav_tab_scaled[rgb] = (float **)calloc(h, sizeof(float *));
	  for(i = 0; i < h; i++) { /* for each freq band */
	       int wl = pitch[i];
	       float scale = (float)TABLE_SIZE / (float)wl;
	       int s;
	       int phase = random() % TABLE_SIZE;
	       wav_tab_scaled[rgb][i] = (float *)calloc(wl,sizeof(float));
	       for(s = 0; s < wl; s++) {
		    int samp = (int)((float)s * scale);
		    wav_tab_scaled[rgb][i][s] = wav_tab[rgb][(samp + phase) % TABLE_SIZE];
	       }
	  }
     }

     /* init some thread attrs */
     pthread_attr_init(&thread_attr);
     pthread_attr_setdetachstate(&thread_attr,PTHREAD_CREATE_JOINABLE);

     /*
       ampScale = 0.0075 / (float)h;
     */
     ampScale /= (float)h;
     /* now produce the sound */
     for(i = 1; i < w-1; i++) { /* for each time segment, */
	  for(j = 0; j < h; j++) { /* for each frequency band */
	       /* compute the energy based on color */
	       int c;
	       uint8 rgb[3];
	       uint32 pixel;
	       uint8 rgb_prev[3];
	       uint32 pixel_prev;

	       pixel = image[(j * w) + i];

	       rgb[0] = (uint8)((pixel & 0x00FF0000) >> 16);
	       rgb[1] = (uint8)((pixel & 0x0000FF00) >> 8);
	       rgb[2] = (uint8)((pixel & 0x000000FF) >> 0);

	       pixel_prev = image[(j * w) + (i - 1)];

	       rgb_prev[0] = (uint8)((pixel_prev & 0x00FF0000) >> 16);
	       rgb_prev[1] = (uint8)((pixel_prev & 0x0000FF00) >> 8);
	       rgb_prev[2] = (uint8)((pixel_prev & 0x000000FF) >> 0);

	       for(c = 0; c < 3; c++) {
		    energy_prev[c][j] = energy[c][j];
		    if(!bischMode) {
			 /* energy for each wave type is inverse color intensity */
			 energy[c][j] = ampScale * freqResponse[j] *
			      (255 - rgb[c]);
		    } else { /* bisch mode */
			 /* energy for each wave type is absolute color difference
			    from previous pixel */
			 energy[c][j] = ampScale * freqResponse[j] *
			      abs((int)rgb[c] - (int)rgb_prev[c]);
		    }
	       }
	  }
	  /* map */
	  for(k = 0; k < nThreads; k++) {
	    ts[k].fi = k * (h / nThreads);
	    ts[k].fn = (k+1) * (h / nThreads);
	    ts[k].pitch = pitch;
	    ts[k].energy = energy;
	    ts[k].energy_prev = energy_prev;
	    ts[k].offset = offset;
	    ts[k].pan = pan;
	    ts[k].ss = cs;
	    ts[k].ns = block;
	    ts[k].buf = &result[k * block * N_CHANNELS];
	    if(k == 0) {
	      do_timeSegment(&ts[k]);
	    } else {
	      pthread_create(&threads[k],&thread_attr,do_timeSegment,(void *)&(ts[k]));
	    }
	  }
	  /* reduce */
	  for(k = 1; k < nThreads; k++) {
	    void *status;
	    pthread_join(threads[k],status);
	    for(s = 0; s < block * N_CHANNELS; s++) {
	      result[s] += result[s+(k*block*N_CHANNELS)];
	    }
	  }
	  output(out, result, block);
	  cs += block;
	  if(!(i % 10)) {
	       show_progress(i,w-1);
	  }
     }
     final_output(out);
     show_progress(100,100);
     fprintf(stderr,"\n");

/*
     result = do_timeSegment(3, pitch, energy, 0, block);
     output(result,block);
     result = do_timeSegment(3, pitch, energy, block, block);
     output(result,block);
*/

     return 1;
}
