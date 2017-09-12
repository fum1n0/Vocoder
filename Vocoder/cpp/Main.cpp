
#include"../hpp/LPC.hpp"

void Main(){

	Console::Open();

	auto lpc = std::make_shared<LPC>();
	lpc->calc_formant(24);
	//lpc->calc_AllFormant();
	
	Console::Close();

}
