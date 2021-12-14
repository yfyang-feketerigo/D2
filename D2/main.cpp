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
#include <boost/filesystem.hpp>

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
	namespace fs = boost::filesystem;
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
		fs::path fs_iffpath_t;
		fs::path fs_iffpath_t_;
		string pairstyle_t;
		string pairstyle_t_;
		string boxtype_t;
		string boxtype_t_;

		po::options_description if_options("infile options");
		if_options.add_options()
			("ifname_t,t", po::value<string>(&ifname_t), "input file name at time t")
			("ifname_t-dt,0", po::value<string>(&ifname_t_), "input file name at time t-dt")
			("ifpath,I", po::value<string>(&ifpath)->default_value("./"), "input file path")
			("pairstyle_t", po::value<string>(&pairstyle_t)->default_value("none"), "input pair style of data file at time t\n allowed args:'none' 'pair' 'single'")
			("pairstyle_t-dt", po::value<string>(&pairstyle_t_)->default_value("none"), "input pair style of data file at time t-dt\n allowed args:'none' 'pair' 'single'")
			("boxtype_t", po::value<string>(&boxtype_t)->default_value("tilt"), "input box type of data file at time t\n allowed args:'orthogonal' 'tilt'")
			("boxtype_t_", po::value<string>(&boxtype_t_)->default_value("tilt"), "input box type of data file at time t-dt\n allowed args:'orthogonal' 'tilt'");

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
		{
			fs_iffpath_t = fs::path(ifpath) / ifname_t;
			cout << "input file name at time t: " << fs_iffpath_t.string() << '\n';
			cout << "input pair style at time t: " << pairstyle_t << '\n';
			cout << "input box type at time t: " << boxtype_t << '\n';
		}
		else
			throw exception("inpute file name at time t was not set.\n --help to show help");

		if (vm.count("ifname_t-dt"))
		{
			fs_iffpath_t_ = fs::path(ifpath) / ifname_t_;
			cout << "input file name at time t-dt: " << fs_iffpath_t_.string() << '\n';
			cout << "input pair style at time t-dt: " << pairstyle_t_ << '\n';
			cout << "input box type at time t-dt: " << boxtype_t_ << '\n';
		}
		else
			throw exception("input file name at time t-dt was not set.\n --help to show help");

		fs::path fs_offpath = fs::path(ofpath) / ofname;
		cout << "output file: " << fs_offpath << '\n';

		cout << '\n';


		using BoxType = Configuration::Configuration::BoxType;
		using PairStyle = Configuration::Configuration::PairStyle;
		BoxType bt_t = Configuration::Configuration::string_to_BoxType(boxtype_t);
		PairStyle ps_t = Configuration::Configuration::string_to_PairStyle(pairstyle_t);
		BoxType bt_t_ = Configuration::Configuration::string_to_BoxType(boxtype_t_);
		PairStyle ps_t_ = Configuration::Configuration::string_to_PairStyle(pairstyle_t_);

		D2::Configuration_neighbours config_t(fs_iffpath_t.string(), rcut, bt_t, ps_t);
		D2::Configuration_neighbours config_t_(fs_iffpath_t_.string(), rcut, bt_t_, ps_t_);

		config_t.sort_particle();
		config_t_.sort_particle();
		config_t.sort_neighbours_as_center_pid();
		config_t.sort_neighbours_as_center_pid();
		//auto& vec_neighbours = config_t.get_neighbours();
		const auto& vec_pa_t_ = config_t_.get_particle();

		vector<double> vecd_D2;
		vecd_D2.resize(vec_pa_t_.size());
		for (size_t i = 0; i < vec_pa_t_.size(); i++)
		{
			vecd_D2[i] = 0.;
			const Particle& pa = vec_pa_t_[i];
			D2::D2 d2_pa(&pa, rcut, &config_t, &config_t_);
			vecd_D2[i] = d2_pa.get_D2();
		}
		mkdir(ofpath);
		config_t_.para_to_dump(fs_offpath.string(), { "D2" }, { vecd_D2 });
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
