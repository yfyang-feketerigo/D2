#include "D2.h"
#include "particle.h"
#include "configuration.h"
#include "input.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>


void mkdir(std::string path)
{
	boost::filesystem::path bpath(path);
	if (!boost::filesystem::exists(path))
	{
		boost::filesystem::create_directories(bpath);
	}
}

int main(int ac, char* av[])
{
	boost::timer::cpu_timer timer;
	namespace po = boost::program_options;
	using namespace std;
	try
	{
		po::options_description help_options("help options");
		help_options.add_options()
			("help", "produce help message");


		double rcut;
		po::options_description rcut_options("rcut");
		rcut_options.add_options()
			("rcut,r", po::value<double>(&rcut), "cutoff distance");

		string ifname_t;
		string ifname_t_;
		string ifpath;
		po::options_description if_options("infile options");
		if_options.add_options()
			("ifname_t,t", po::value<string>(&ifname_t), "input file name at time t")
			("ifname_t-dt,T", po::value<string>(&ifname_t_), "input file name at time t-dt")
			("ifpath,I", po::value<string>(&ifpath)->default_value("./"), "input file path");
		po::options_description of_options("output file options");

		string ofname;
		string ofpath;
		of_options.add_options()
			("ofname,o", po::value<string>(&ofname)->default_value("D2"), "D2 output file name")
			("ofpath,O", po::value<string>(&ofpath)->default_value("./"), "D2 output file path");

		po::options_description allowed_options("Allowed options");
		allowed_options.add(help_options).add(rcut_options).add(if_options).add(of_options);

		po::variables_map vm;
		po::store(po::parse_command_line(ac, av, allowed_options), vm);
		po::notify(vm);
		if (vm.count("help"))
		{
			cout << allowed_options << "\n";
			return 1;
		}

		if (vm.count("rcut"))
			cout << "rcut: " << rcut << '\n';
		else
			throw exception("rcut was not set!");

		if (vm.count("ifname_t"))
			cout << "input file name at time t: " << vm["ifname_t"].as<string>() << '\n';
		else
			throw exception("inpute file name at time t was not set.\n --help to show help");

		if (vm.count("ifname_t-dt"))
			cout << "input file name at time t-dt: " << vm["ifname_t-dt"].as<string>() << '\n';
		else
			throw exception("input file name at time t-dt was not set.\n --help to show help");

		cout << "input file path: " << vm["ifpath"].as<string>() << '\n';
		cout << "output file name: " << vm["ofname"].as<string>() << '\n';
		cout << "output file path: " << vm["ofpath"].as<string>() << '\n';
		cout << '\n';


		using BoxType = Configuration::Configuration::BoxType;
		using PairStyle = Configuration::Configuration::PairStyle;
		//typedef Configuration::Configuration::BoxType BoxType;
		//typedef Configuration::Configuration::PairStyle PairStyle;

		string fname_t = ifpath + '/' + ifname_t;
		string fname_t_ = ifpath + '/' + ifname_t_;

		D2::Configuration_neighbours config_t(fname_t, rcut, BoxType::tilt, PairStyle::none);
		D2::Configuration_neighbours config_t_(fname_t_, rcut, BoxType::tilt, PairStyle::none);

		config_t.sort_particle();
		config_t_.sort_particle();
		config_t.sort_neighbours_as_center_pid();
		config_t.sort_neighbours_as_center_pid();
		auto& vec_neighbours = config_t.get_neighbours();
		const auto& vec_pa_t = config_t.get_particle();

		vector<double> vecd_D2;
		vecd_D2.resize(vec_pa_t.size());
		for (size_t i = 0; i < vec_pa_t.size(); i++)
		{
			vecd_D2[i] = 0.;
			const Particle& pa = config_t.get_particle()[i];
			D2::D2 d2_pa(&pa, rcut, &config_t, &config_t_);
			vecd_D2[i] = d2_pa.get_D2();
		}
		mkdir(ofpath);
		string offname = ofpath + '/' + ofname;
		config_t.para_to_dump(offname, { "D2" }, { vecd_D2 });
		cout << "done! time info: " << timer.format();
		cout << '\n';
	}
	catch (const std::exception& e)
	{
		cerr << "error: " << e.what() << '\n';
	}
	catch (...) {
		cerr << "Exception of unknown type!\n";
	}
	return 0;
}
