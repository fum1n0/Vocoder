
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
	// 0.0ˆÈã1.0–¢–‚Ì’l‚ğ“™Šm—¦‚Å”­¶‚³‚¹‚é
	std::uniform_real_distribution<> dist(-1.0, 1.0);


	/*Wave white(Frame_L);
	for (int i = 0; i < Frame_L; i++)white[i] = Waving::DoubleToSample(dist(engine));
	white.saveWAVE(L"whitenoise.wav");*/


	for (int i = 0; i < Frame_L; i++)whitenoise[i] = dist(engine);
	//for (int i = 1; i < Frame_L; i++)whitenoise[i] = (whitenoise[i] + whitenoise[i - 1]) / 2.0;
	//for (int i = 0; i < Frame_L; i++)whitenoise[i] *= sqrt(3.0)*e_rms;

	
	/*Wave white(Frame_L);
	for (int i = 0; i < Frame_L; i++)white[i] = Waving::DoubleToSample(whitenoise[i]);
	white.saveWAVE(L"whitenoise.wav");*/

}


void Whisper::createWhisper() {

	double whisper_amp;
	whisperVoice = std::vector<double>(wav_pre.lengthSample, 0);

	for (long i = 0; i < leg; i++) {

		createWhiteNoise();
		LPC::calc_formant(i);
	
		for (int j = 0; j < Order; j++) {
			//whisperVoice.push_back(whitenoise[j]);
			//whisperVoice[i*Frame_T + j] += whitenoise[j];
			
			whisperVoice[i*Frame_T + j] += wav_pre[i*Frame_T + j].left/32768.0;
		}

		for (int j = Order; j < Frame_L; j++) {
			whisper_amp = 0;
			for (int k = 1; k < a.size(); k++) {
				//whisper_amp -= a[k] * (whitenoise[j - k]);
				//whisper_amp -= a[k] * whisperVoice[i*Frame_T + j-k];

				whisper_amp -= a[k] * wav_pre[i*Frame_T + j - k].left / 32768.0;
			}

			whisperVoice[i*Frame_T + j] += whisper_amp;
		}

	}

	calc_Normalization(whisperVoice);

	voice = Wave(whisperVoice.size());
	for (long i = 0; i < whisperVoice.size(); i++)voice[i] = Waving::DoubleToSample(whisperVoice[i]);

	voice.saveWAVE(Widen(whisperName + ".wav"));

}