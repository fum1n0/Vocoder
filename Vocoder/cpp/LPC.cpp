#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include"../hpp/LPC.hpp"


LPC::LPC() {
	setLPC(1024, 128);
}


LPC::LPC(int frameL, int frameT) {

	if ((frameL & (frameL - 1)) != 0 || (frameT & (frameT - 1)) != 0) {
		std::cout << "Constructor Error" << std::endl;
		std::cout << "Set Default Number" << std::endl;
		this->LPC::LPC();
	}else setLPC(frameL, frameT);
	
}


void LPC::setLPC(int frameL, int frameT) {
		
	while (sound.isEmpty() ) {
		Optional<String> state = Dialog::GetOpenSound();
		if (!state.has_value()) {
			std::cout << "Sound Select Error" << std::endl;
			continue;
		}
		sound = Sound(state.value());
		Optional<AudioProperty> audio = Audio::GetProperty(state.value());
		name = audio->title.narrow();
		name = name.substr(0, name.size() - 4);
	}

	fileFormantFreq_Band = name + "_Formant.txt";
	writing_FormantFreq_Band.open(fileFormantFreq_Band, std::ios::out);

	Frame_L = frameL;
	Frame_T = frameT;

	wav_zero = Wave(1); // 0-pound

	wav_base = sound.getWave();
	wav_base = erase_zeroAmp(wav_base);
	length = wav_base.lengthSample;
	leg = (long)(length / Frame_T);

	wav_zero[0] = Waving::DoubleToSample(0.0);

	for (int m = 0; m < Frame_L; m++)wav_base.push_back(wav_zero[0]);

	wav_pre = wav_base;
	wav_ana = Wave(Frame_L);
	
	// 高音強調
	for (size_t m = 1; m < length; m++)
		wav_pre[m] = Waving::DoubleToSample(((double)wav_base[m].left / 32768.0) - ((double)wav_base[m - 1].left / 32768.0) * pre_emphasis);



	// Order = (int)(sound.samplingRate() / 1000) + 2;

	signal = std::vector<double>(Frame_L, 0);

	re = std::vector<double>(Frame_L, 0);;
	im = std::vector<double>(Frame_L, 0);
	
	amp = std::vector<double>(Frame_L, 0); // Amp
	power = std::vector<double>(Frame_L, 0);//power
	dbv = std::vector<double>(Frame_L, 0); // dBV

	r = std::vector<double>(Frame_L, 0); // 自己相関関数

	a = std::vector<double>(Order + 1, 0); // 線形予測係数
	b = std::vector<double>(Order + 1, 0); // 誤差信号係数
	
	e = std::vector<double>(Frame_L, 0); // 残差信号
	res_auto = std::vector<double>(Frame_L, 0); // 残差信号の自己相関関数

	lpc_gain = std::vector<double>(Frame_L, 0);
	lpc_dbv = std::vector<double>(Frame_L, 0);;

	formant_freq = std::vector<std::vector<double>>(leg + 1, std::vector<double>());
	

}

LPC::~LPC() {
	writing_FormantFreq_Band.close();
}

Wave LPC::erase_zeroAmp(Wave& wav) {

	Wave erase;
	bool flag = true;
	for (size_t i = 1; i < wav.lengthSample - 1; i++) {
		flag = true;

		if (wav[i - 1].left == wav[i].left) {
			if (wav[i + 1].right == wav[i + 1].right) {
				flag = false;
			}
		}

		if (flag)erase.push_back(wav[i]);

	}

	erase.saveWAVE(L"erase.wav");

	return erase;


}

void LPC::calc_AllFormant() {

	for (int i = 0; i < leg; i++)calc_formant(i);

}


void LPC::calc_formant(int indexFrame) {

	nowIndexFrame = indexFrame;

	LPC::init();
	LPC::hanning_execute(indexFrame);
	LPC::preEmphasis();
	LPC::calc_ACF_FFT();
	LPC::calc_Levinson_Durbin();

	LPC::calc_error();
	//LPC::calc_PitchFreq();
	
	LPC::calc_lpc_gain();
	//LPC::calc_Normalization(lpc_gain);
	LPC::calc_lpc_dBV();
	
	//LPC::calc_roots();
	//LPC::count_formant(indexFrame);

}


void LPC::init() {
	signal = std::vector<double>(Frame_L, 0);
	y = std::vector<double>();
	re = std::vector<double>(Frame_L, 0);
	im = std::vector<double>(Frame_L, 0);
	amp = std::vector<double>(Frame_L, 0);
	power = std::vector<double>(Frame_L, 0);
	dbv = std::vector<double>(Frame_L, 0);
	
	lpc_gain = std::vector<double>(Frame_L, 0);
	lpc_dbv = std::vector<double>(Frame_L, 0);

	a = std::vector<double>(Order + 1, 0);
	b = std::vector<double>(Order + 1, 0);
	
	r = std::vector<double>(Frame_L, 0);
	e = std::vector<double>(Frame_L, 0);
	res_auto = std::vector<double>(Frame_L, 0);
	e_rms = 0.0;

	wav_ana = Wave(Frame_L);
}

void LPC::hanning_execute(int indexFrame) {

	for (int k = 0; k < Frame_L; k++) {
		wav_ana[k] = wav_base[indexFrame*Frame_T + k];
		signal[k] = (wav_base[indexFrame*Frame_T + k].left/32768.0 ) * (0.5 - 0.5*cos(2.0 * Pi * k / (double)Frame_L)); // ハニング補正無し
	}
	

}


void LPC::preEmphasis() {

	std::string file_signal = "signal.txt";
	std::ofstream writing_signal;
	writing_signal.open(file_signal, std::ios::out);
	
	for (int k = 1; k < Frame_L; k++) {
		signal[k] = signal[k]-pre_emphasis*signal[k-1];
	}
	
	for (int i = 0; i < Frame_L;i++)writing_signal << i << " " << signal[i] << std::endl;

	writing_signal.close();
	
}


void LPC::fft_excute(std::vector<double>& signal_ptr, std::vector<double>& re_ptr, std::vector<double>& im_ptr, int inv) {

	std::vector<double>sig = signal_ptr;

	double w;
	std::vector<double>wr;
	std::vector<double>wi;
	w = (double)inv * 2.0 * Pi / (double)sig.size();
	wr = std::vector<double>((int)sig.size(), 0);
	wi = std::vector<double>((int)sig.size(), 0);
	for (int i = 0; i < (int)sig.size(); i++) {
		wr[i] = cos(w * i);
		wi[i] = -sin(w*i);
	}


	long i = 0;

	if (inv > 0) {
		for (int j = 1; j < (int)sig.size() - 1; j++) {
			for (int k = (int)sig.size() >> 1; k > (i ^= k); k >>= 1); // k=2^{n-1}; k > i = i xor k; k/2
			if (j < i) {
				std::swap(sig[i], sig[j]);
			}
		}
		re_ptr = sig;

	}
	else {
		for (int j = 1; j < (int)sig.size() - 1; j++) {
			for (int k = (int)sig.size() >> 1; k > (i ^= k); k >>= 1); // k=2^{n-1}; k > i = i xor k; k/2
			if (j < i) {
				std::swap(re_ptr[i], re_ptr[j]);
				std::swap(im_ptr[i], im_ptr[j]);
			}
		}
	}


	double xr1, xr2, xi1, xi2;
	int s, t, turn1, turn2, turn3;

	for (int rep = 2; rep <= (int)sig.size(); rep *= 2) { // 一番外側のrep, 2DFT -> 4DFT -> ... -> NDFT
		s = (int)sig.size() / rep; // ブロック数
		t = rep / 2; // 1ブロックが保有する数の半分
		for (int u = 0; u < s; u++) { // NDFTの中のブロック数
			for (int v = 0; v < t; v++) { // ブロックの半分まで

				turn1 = v + rep*u;
				turn2 = s*v;
				turn3 = s*(v + t);

				xr1 = re_ptr[turn1];
				xr2 = re_ptr[turn1 + t];
				xi1 = im_ptr[turn1];
				xi2 = im_ptr[turn1 + t];


				re_ptr[turn1] = xr1 + xr2 * wr[turn2] - xi2 * wi[turn2];
				im_ptr[turn1] = xi1 + xi2 * wr[turn2] + xr2 * wi[turn2];

				re_ptr[turn1 + t] = xr1 + xr2 * wr[turn3] - xi2 * wi[turn3];
				im_ptr[turn1 + t] = xi1 + xi2 * wr[turn3] + xr2 * wi[turn3];


			}
		}
	}

	if (inv > 0) { // フーリエ変換
		for (int k = 0; k < (int)re_ptr.size(); k++) {
			re_ptr[k] /= (double)re_ptr.size();
			im_ptr[k] /= (double)im_ptr.size();
		}
	}
	else signal_ptr = re_ptr; // 逆フーリエ変換


}


void LPC::calc_ACF_FFT() {

	std::vector<double> sig = signal;
	std::vector<double> sigHanning = signal;

	for (auto& i : sigHanning)i = i*2.0;

	fft_excute(sigHanning, re, im, 1); // フーリエ変換

	std::string file_dbv = "signal_dbv.txt";
	std::ofstream writing_dbv;
	writing_dbv.open(file_dbv, std::ios::out);
	for (int i = 0; i < Frame_L; i++) {
		amp[i] = sqrt(re[i] * re[i] + im[i] * im[i]);
		dbv[i] = 20 * log10(amp[i]);
		writing_dbv << (double)sound.samplingRate()*i / (double)Frame_L << " " << dbv[i] << std::endl;
	}
	writing_dbv.close();

	for (int i = 0; i < Frame_L; i++)sig.push_back(0.0); // 0-padding

	re = std::vector<double>((int)sig.size(), 0);
	im = std::vector<double>((int)sig.size(), 0);
	power = std::vector<double>((int)sig.size(), 0);

	fft_excute(sig, re, im, 1); // フーリエ変換

	for (int i = 0; i < (int)power.size(); i++)power[i] = (re[i] * re[i] + im[i] * im[i]);
	
	re = power;
	im = std::vector<double>((int)re.size(), 0);
	fft_excute(sig, re, im, -1); // 逆フーリエ変換

	for (int i = 0; i < (int)r.size(); i++)r[i] = sig[i];

}

void LPC::calc_ACF_FFT(std::vector<double>& r_ptr, std::vector<double>& signal_ptr) {
	std::vector<double> sig = signal_ptr;

	for (int i = 0; i < Frame_L; i++)sig.push_back(0.0); // 0-padding

	re = std::vector<double>((int)sig.size(), 0);
	im = std::vector<double>((int)sig.size(), 0);
	
	fft_excute(sig, re, im, 1); // フーリエ変換
		
	for (int i = 0; i < (int)power.size(); i++)power[i] = (re[i] * re[i] + im[i] * im[i]);
	
	re = power;
	im = std::vector<double>((int)re.size(), 0);
	fft_excute(sig, re, im, -1); // 逆フーリエ変換

	for (int i = 0; i < (int)r_ptr.size(); i++)r_ptr[i] = sig[i];
}


void LPC::calc_Levinson_Durbin() {
	int k, j;

	a[0] = b[0] = 1.0;
	a[1] = -r[1] / r[0];
	b[1] = r[0] + r[1] * a[1];
	double lambda = -r[1] / r[0];

	std::vector<double>U(a.size() + 1);
	std::vector<double>V(a.size() + 1);


	for (k = 1; k < Order; k++) {
		lambda = 0.0;
		for (j = 0; j < k + 1; j++) {
			lambda -= a[j] * r[k + 1 - j];
		}
		lambda /= b[k];

		U[0] = 1.0;
		for (j = 1; j < k + 1; j++)U[j] = a[j];
		U[j] = 0.0;

		V[0] = 0.0;
		for (j = k; j > 0; j--)V[k - j + 1] = a[j];
		V[k - j + 1] = 1.0;
		
		for (int s = 0; s <= k+1; s++) {
			a[s] = U[s] + V[s] * lambda;
		}

		b[k + 1] = b[k] * (1.0 - lambda * lambda);
	}
}


void LPC::calc_error() {
	// 線形予測
	double tmp;
	// 補正あり,正規化
		
	for (int k = 0; k < Frame_L; k++) {
		tmp = 0;
		if (k < Order)y.push_back(signal[k]);
		
		else {
			for (int j = 1; j < (int)a.size(); j++) {
				tmp -= a[j] * (signal[k - j]);
			}
			y.push_back(tmp);
		}
	}

	//残差信号
	for (int k = 0; k < Frame_L; k++) {
			
		
		if (Order < k) {
			tmp = 0.0;
			for (int j = 0; j < (int)a.size(); j++) {
				tmp += a[j] * (signal[k - j]);
			}
			//e[k] = tmp;
			e[k] = tmp * (0.5 - 0.5*cos(2 * Pi*k / Frame_L)); // ハニング窓
		}
		
		e_rms += e[k]*e[k];
	}

	e_rms = sqrt(e_rms/(double)Frame_L);


	// output file
	std::string file_y = "y.txt";
	std::ofstream writing_y;
	std::string file_e = "e.txt";
	std::ofstream writing_e;	
		
	writing_y.open(file_y, std::ios::out);
	writing_e.open(file_e, std::ios::out);

	for (int i = 0; i < Frame_L; i++) {
		writing_y << i << " " << y[i] << std::endl;
		writing_e << i << " " << e[i] << std::endl;
	}

	
	writing_y.close();
	writing_e.close();

}


void LPC::calc_lpc_gain() {

	double re_tmp, im_tmp;

	std::ofstream writingGain;
	writingGain.open("lpc_gain.txt", std::ios::out);

	for (int i = 0; i < (int)lpc_gain.size(); i++) {
		re_tmp = im_tmp = 0;

		for (int j = 0; j < (int)a.size(); j++) {

			re_tmp +=  (a[j] * cos(2 * Pi*i*j / (double)Frame_L));
			im_tmp +=  (a[j] * (-1) * sin(2 * Pi*i*j / (double)Frame_L));

		}
		
		lpc_gain[i] = 1.0 / (sqrt(re_tmp*re_tmp + im_tmp*im_tmp));
		writingGain << lpc_gain[i] << std::endl;
	}
	writingGain.close();

}


void LPC::calc_Normalization(std::vector<double>& num) {

	double max = *std::max_element(num.begin(), num.end());
	for (int i = 0; i < (int)num.size(); i++)num[i] /= max;

}


void LPC::calc_lpc_dBV() {

	std::string file_dbv = "lpc_dbv.txt";
	std::ofstream writing_dbv;
	writing_dbv.open(file_dbv, std::ios::out);
	
	// calc lpc dBV
	for (int i = 0; i < (int)lpc_gain.size(); i++)lpc_dbv[i] = 20 * log10(lpc_gain[i] * e_rms);
	
	// move lpc dBV graph
	/*double lpc_dbv_dif = *std::max_element(lpc_dbv.begin(), lpc_dbv.end()) - *std::max_element(dbv.begin(), dbv.end());;
	for (int i = 0; i < (int)lpc_gain.size(); i++) {
		writing_dbv << (double)sound.samplingRate()*i / (double)Frame_L << " " << lpc_dbv[i] - lpc_dbv_dif << std::endl;
	}*/

	for (int i = 0; i < (int)lpc_gain.size(); i++) {
		writing_dbv << (double)sound.samplingRate()*i / (double)Frame_L << " " << lpc_dbv[i] << std::endl;
	}


	writing_dbv.close();

}


void LPC::count_formant(int indexFrame) {
		
	int frek;

	for (int l = 0; l < (int)freq_band.size(); l++) {
		//writing_FormantFreq_Band << freq_band[l].first << " " << freq_band[l].second << std::endl;
		frek = (int)(freq_band[l].first / 1000);

		if (band_min <= freq_band[l].second && freq_band[l].second <= band_max ) {
			formant_freq[indexFrame].push_back(freq_band[l].first);
			writing_FormantFreq_Band << freq_band[l].first << " " << freq_band[l].second << std::endl;
		}

	}
	writing_FormantFreq_Band << std::endl;

}


void LPC::calc_roots() {

	std::vector<std::vector<double>>companion(a.size() - 1, std::vector<double>(a.size() - 1, 0));
	reverse(a.begin(), a.end());
	for (int i = 0; i < (int)companion.size(); i++) {
		if (0 < i)companion[i][i - 1] = 1;
		companion[i][companion[i].size() - 1] = -a[i];
	}


	Eigen::MatrixXf eigen_companion = eigen_matrix(companion);
	//std::cout << eigen_companion << std::endl;

	Eigen::MatrixXd eigen_companion_d(companion.size(), companion[0].size());
	for (int i = 0; i < companion.size(); ++i)
		eigen_companion_d.col(i) = Eigen::VectorXd::Map(&companion[i][0], companion[0].size());
	//std::cout << eigen_companion_d << std::endl;


	Eigen::Matrix<std::complex<double>, 46, 46> plex;

	for (int i = 0; i < 46; i++) {
		for (int j = 0; j < 46; j++) {
			plex(j, i) = std::complex<double>(eigen_companion(j, i), 0); // jが行,iが列
		}
	}


	Eigen::ComplexEigenSolver<Eigen::Matrix<std::complex<double>, 46, 46>>s(plex);

	Eigen::Matrix<std::complex<double>, 46, 1> ve = s.eigenvalues();

	freq_band = std::vector<std::pair<double, double>>();

	double img ;
	double real;
	double freq;

	for (int i = 0; i < ve.size(); i++) {
		img = ve(i, 0).imag();
		real = ve(i, 0).real();
		if (img > 0) {
			freq = 0;

			if (real > 0)freq = (double)sound.samplingRate()*(atan(img / real)) / (2 * Pi);
			else freq = (double)sound.samplingRate()*(Pi + atan(img / real)) / (2 * Pi);

			double ban = -1 * (double)sound.samplingRate()*log(abs(sqrt(img*img + real*real))) / Pi;
			freq_band.push_back(std::pair<double, double>(freq, ban));

		}
	}

	sort(freq_band.begin(), freq_band.end());

}


void LPC::calc_PitchFreq() {

	LPC::calc_error();
	LPC::calc_ACF_FFT(res_auto, e);
	
	std::vector<double>::iterator maxIt = std::max_element(res_auto.begin()+1, res_auto.end());
	size_t maxIndex = std::distance(res_auto.begin(), maxIt);

	pitch_freq = sound.samplingRate() / (double)maxIndex;
	std::cout <<"Pitch Freq: "<<pitch_freq <<"[Hz]"<< std::endl;

}