
#include"../hpp/Vocoder.hpp"

Vocoder::Vocoder() :LPC(Frame_L, Frame_L) {
	setWhisperParameter();
}


Vocoder::Vocoder(int frameL, int frameT) : LPC(frameL, frameT) {
	setWhisperParameter();
}

void Vocoder::setWhisperParameter() {
	whisperName = name + "_whisper";
	robotName = name + "_robot";
	whitenoise = std::vector<double>(Frame_L, 0);
	sawtooth = std::vector<double>(Frame_L, 0);
	time = 1.0/(double)sound.samplingRate();
}


void Vocoder::createWhiteNoise() {

	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());
	std::uniform_real_distribution<> dist(-1.0, 1.0);
		
	for (int i = 0; i < Frame_L; i++)whitenoise[i] = dist(engine)*e_rms;
		
}


void Vocoder::createSawtooth() {
	for (int i = 0; i < Frame_L; i++)
		sawtooth[i] = 2 * (260.0*i*time - floor(260.0*i*time + 0.5))*e_rms;

}



void Vocoder::createWhisper() {

	whisperVoice = std::vector<double>(wav_pre.lengthSample, 0);

	std::vector<double>amp_e(wav_pre.lengthSample, 0);

	for (long i = 0; i < leg; i++) {

		
		LPC::calc_formant(i);
		createWhiteNoise();

		for (long j = 0; j<Frame_L; j++) {
			for (long k = 1; k<a.size(); k++) {
				if (j >= k) whitenoise[j] -= a[k] * whitenoise[j - k];
			}
			whisperVoice[i*Frame_T + j] += whitenoise[j];
			amp_e[i*Frame_T + j] += e[j];
		}
	}

	for (long i = 1; i < whisperVoice.size(); i++)whisperVoice[i] += de_emphasis*whisperVoice[i - 1];
	//calc_Normalization(whisperVoice);
	
	voice = Wave(whisperVoice.size());
	for (long i = 0; i < whisperVoice.size(); i++)voice[i] = Waving::DoubleToSample(whisperVoice[i]);
	voice.saveWAVE(Widen(whisperName + ".wav"));
		
	for (long i = 1; i < amp_e.size(); i++)amp_e[i] = (amp_e[i] + amp_e[i - 1]) / 2.0;
	//calc_Normalization(amp_e);
	Wave error(amp_e.size());
	for(long i=0;i<amp_e.size();i++)error[i] = Waving::DoubleToSample(amp_e[i]);
	error.saveWAVE(Widen(whisperName + "_e.wav"));

}


void Vocoder::createRobot() {
	robotVoice = std::vector<double>(wav_pre.lengthSample, 0);

	std::vector<double>amp_e(wav_pre.lengthSample, 0);

	for (long i = 0; i < leg; i++) {

		
		LPC::calc_formant(i);
		createSawtooth();

		for (long j = 0; j<Frame_L; j++) {
			for (long k = 1; k<a.size(); k++) {
				if (j >= k) sawtooth[j] -= a[k] * sawtooth[j - k];
			}
			robotVoice[i*Frame_T + j] += sawtooth[j];
			amp_e[i*Frame_T + j] += e[j];
		}
	}

	for (long i = 1; i < robotVoice.size(); i++)robotVoice[i] += de_emphasis*robotVoice[i - 1];
	calc_Normalization(robotVoice);

	voice = Wave(robotVoice.size());
	for (long i = 0; i < robotVoice.size(); i++)voice[i] = Waving::DoubleToSample(robotVoice[i]);
	voice.saveWAVE(Widen(robotName + ".wav"));

	for (long i = 1; i < amp_e.size(); i++)amp_e[i] = (amp_e[i] + amp_e[i - 1]) / 2.0;
	calc_Normalization(amp_e);
	Wave error(amp_e.size());
	for (long i = 0; i<amp_e.size(); i++)error[i] = Waving::DoubleToSample(amp_e[i]);
	error.saveWAVE(Widen(robotName + "_e.wav"));
}