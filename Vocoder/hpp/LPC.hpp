#pragma once


#include"matrix_util.hpp"


class LPC{

protected:
		
	Sound sound; // �����t�@�C��
	std::string name; // �����t�@�C����

	Wave wav_base; // ���g�`
	Wave wav_pre; // ��������������g�`
	Wave wav_zero; // 0-pound

	unsigned long long int length ; // �����T���v����
	int Frame_L = 1024; // �t���[����
	int Frame_T = 128; // �t���[������
	long leg; // ������ = �����T���v����/�t���[������
	long nowIndexFrame; // ����͒��̎���

	double pre_emphasis = 0.92; // ������������
	double de_emphasis = 0.30; // �ቹ��������

	double e_rms = 0.0; // ��敽�ϕ�����

	double pitch_freq; // �s�b�`���g��

	Wave wav_ana; // ��͗p�t���[���������؂�o������
	
	const int Order = 46; // ��͎���
	const unsigned int band_min = 0; // �t�H���}���g�ш扺��
	const unsigned int band_max = 350; // �t�H���}���g�ш���

	std::vector<double>signal; // �t���[�����؂�o���V�O�i��

	std::vector<double>re; // ��������
	std::vector<double>im; // ��������

	std::vector<double>amp; // Amp
	std::vector<double>power; //power
	std::vector<double>dbv; // dBV

	std::vector<double>r; // ���ȑ��֊֐�

	std::vector<double>a; // ���`�\���W��
	std::vector<double>b; // �덷�M���W��
	
	std::vector<double>y; // �\���M��

	std::vector<double>e; // �c���M��
	std::vector<double>res_auto; // �c���M���̎��ȑ��֊֐�

	std::vector<double>lpc_gain; // lpc�̎��g������
	std::vector<double>lpc_dbv; // lpc���g��������dBV

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

	void init();
	void calc_formant(int);
	
	void hanning_execute(int);
	void fft_excute(std::vector<double>&, std::vector<double>&, std::vector<double>&, int);
	void calc_ACF_FFT();
	void calc_ACF_FFT(std::vector<double>&, std::vector<double>&);
	void calc_Levinson_Durbin();
	
	void calc_PitchFreq();
	void calc_error();
	


	void calc_lpc_gain();
	void calc_lpc_dBV();
	void calc_Normalization(std::vector<double>&);
	
	void count_formant(int);
	void calc_roots();
	void calc_AllFormant();

};