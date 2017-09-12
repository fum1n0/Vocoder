#pragma once

#include"LPC.hpp"

class Whisper : public LPC{

private:

	Wave voice;
	std::vector<double>whisperVoice;
	std::vector<double>whitenoise;
	std::string whisperName;

public:

	Whisper();
	Whisper(int, int);

	void setWhisperParameter();
	void createWhisper();
	void createWhiteNoise();

};