

#include"../hpp/Whisper.hpp"

void Main(){

	Console::Open();
	
	auto whisper = std::make_shared<Whisper>(1024,128);
	whisper->createWhisper();
	//whisper->calc_AllFormant();
	//whisper->calc_formant(24);

	Console::Close();

}
