

#include"../hpp/Vocoder.hpp"

void Main(){

	Console::Open();
	
	auto vocoder = std::make_shared<Vocoder>(1024,512);
	//vocoderr->createWhisper();
	vocoder->createRobot();
	Console::Close();

}
