#pragma once


#include"matrix_util.hpp"


class LPC{

private:
		
	Sound sound; // 音声ファイル
	std::string name; // 音声ファイル名

	Wave wav_base; // 元波形
	Wave wav_pre; // 高音強調処理後波形
	Wave wav_zero; // 0-pound

	unsigned long long int length ; // 音声サンプル数
	int Frame_L; // フレーム長
	int Frame_T; // フレーム周期
	long leg; // 周期数 = 音声サンプル数/フレーム周期
	long nowIndexFrame; // 現解析中の周期

	Wave wav_ana;
	
	const int Order = 46; // 解析次数
	const unsigned int band_min = 0;
	const unsigned int band_max = 350;

	std::vector<double>signal; // フレーム長切り出しシグナル

	std::vector<double>re; // 実数部分
	std::vector<double>im; // 虚数部分
	std::vector<double>freq; // 

	std::vector<double>amp; // Amp
	std::vector<double>power; //power
	std::vector<double>dbv; // dBV

	std::vector<double>r; // 自己相関関数

	std::vector<double>a; // 線形予測係数
	std::vector<double>b; // 誤差信号係数
	std::vector<double>root; // 根

	std::vector<double>y; // 予測信号

	std::vector<double>e; // 残差信号
	std::vector<double>res_auto; // 残差信号の自己相関関数

	std::vector<double>lpc_gain; // lpcの周波数応答
	std::vector<double>lpc_dbv; // lpc周波数応答のdBV

	std::vector<std::vector<double>> formant_freq;
	std::vector<std::pair<double, double>>freq_band;
		
	std::string fileFormantFreq_Band;
	std::ofstream writing_FormantFreq_Band;
	

public:

	LPC();
	LPC(int, int);
	~LPC();

	void setLPC(int, int);

	Wave erase_zeroAmp(Wave&);

	void calc_formant(int);
	void init();

	void hanning_execute(int);
	void fft_excute(std::vector<double>&, std::vector<double>&, std::vector<double>&, int);
	void calc_ACF_FFT();
	void calc_ACF_FFT(std::vector<double>&, std::vector<double>&);
	void calc_Levinson_Durbin();
	void calc_error();
	void calc_lpc_gain();
	void calc_lpc_dBV();
	void calc_Normalization(std::vector<double>&);
	void count_formant(int);
	void calc_roots();

	void calc_AllFormant();

};