

#include"../hpp/Whisper.hpp"

void Main(){

	Console::Open();
	
	auto whisper = std::make_shared<Whisper>(1024,512);
	whisper->createWhisper();
	
	Console::Close();

}
