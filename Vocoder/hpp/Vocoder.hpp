#pragma once

#include"LPC.hpp"

class Vocoder : public LPC{

private:

	Wave voice;
	std::vector<double>whisperVoice;
	std::vector<double>robotVoice;
	std::vector<double>whitenoise;
	std::vector<double>sawtooth;
	std::string whisperName;
	std::string robotName;
	double time;

public:

	Vocoder();
	Vocoder(int, int);

	void setWhisperParameter();
	void createWhisper();
	void createWhiteNoise();
	void createSawtooth(); 
	void createRobot();
};