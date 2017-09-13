
#include"../hpp/Whisper.hpp"

Whisper::Whisper() :LPC(Frame_L, Frame_L) {
	setWhisperParameter();
}


Whisper::Whisper(int frameL, int frameT) : LPC(frameL, frameT) {
	setWhisperParameter();
}

void Whisper::setWhisperParameter() {
	whisperName = name + "_whisper";
	whitenoise = std::vector<double>(Frame_L, 0);
}


void Whisper::createWhiteNoise() {

	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());
	std::uniform_real_distribution<> dist(-1.0, 1.0);
		
	for (int i = 0; i < Frame_L; i++)whitenoise[i] = dist(engine)*e_rms;
		
}


void Whisper::createWhisper() {

	whisperVoice = std::vector<double>(wav_pre.lengthSample, 0);

	std::vector<double>amp_e(wav_pre.lengthSample, 0);

	for (long i = 0; i < leg; i++) {

		createWhiteNoise();
		LPC::calc_formant(i);
		
		for (long j = 0; j<Frame_L; j++) {
			for (long k = 1; k<a.size(); k++) {
				if (j >= k) whitenoise[j] -= a[k] * whitenoise[j - k];
			}
			whisperVoice[i*Frame_T + j] += whitenoise[j];
			amp_e[i*Frame_T + j] += e[j];
		}
	}

	for (long i = 1; i < whisperVoice.size(); i++)whisperVoice[i] += de_emphasis*whisperVoice[i - 1];
	calc_Normalization(whisperVoice);
	
	voice = Wave(whisperVoice.size());
	for (long i = 0; i < whisperVoice.size(); i++)voice[i] = Waving::DoubleToSample(whisperVoice[i]);
	voice.saveWAVE(Widen(whisperName + ".wav"));
		
	for (long i = 1; i < amp_e.size(); i++)amp_e[i] = (amp_e[i] + amp_e[i - 1]) / 2.0;
	calc_Normalization(amp_e);
	Wave error(amp_e.size());
	for(long i=0;i<amp_e.size();i++)error[i] = Waving::DoubleToSample(amp_e[i]);
	error.saveWAVE(Widen(whisperName + "_e.wav"));

}