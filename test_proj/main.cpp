#include "../D2/D2.h"
#include "../D2/particle.h"
#include "../D2/configuration.h"
#include "../D2/input.h"
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

int main()
{
	using namespace std;
	try
	{
		using Boxtype = Configuration::Configuration::BoxType;
		using Pairstyle = Configuration::Configuration::PairStyle;
		string testfname = "data.restart.DPD.shear.wi100.5280";
		D2::Configuration_neighbours ptest(testfname, 2.5, Boxtype::tilt, Pairstyle::none);
		string testfname2 = "data.restart.DPD.shear.wi100.6600";
		ptest.sort_particle();
		D2::Configuration_neighbours new_test(ptest);
		D2::Configuration_neighbours ptest2(testfname2, 2.5, Boxtype::tilt, Pairstyle::none);

		cout << ptest.get_neighbours(21).vec_neighbours.size() << '\n';
		cout << new_test.get_neighbours(21).vec_neighbours.size() << '\n';
		ptest = ptest2;
		new_test = ptest2;
		cout << ptest.get_neighbours(21).vec_neighbours.size() << '\n';
		cout << new_test.get_neighbours(21).vec_neighbours.size() << '\n';
	}
	catch (const std::exception& e)
	{
		cerr << "error: " << e.what() << '\n';
		return -1;
	}
	catch (...) {
		cerr << "Exception of unknown type!\n";
		return -10;
	}
	return 0;
}
