#include "D2.h"
#include "particle.h"
#include "configuration.h"
#include "input.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <set>
int main()
{
	using namespace std;
	typedef Configuration::Configuration::BoxType BoxType;
	typedef Configuration::Configuration::PairStyle PairStyle;

	string fname_t = "test_datafile/data.restart.DPD.shear.wi100.6600";
	string fname_t_ = "test_datafile/data.restart.DPD.shear.wi100.5280";
	double r_cut = 2.5;

	D2::Configuration_neighbours config_t(fname_t, r_cut, BoxType::tilt, PairStyle::none);
	D2::Configuration_neighbours config_t_(fname_t_, r_cut, BoxType::tilt, PairStyle::none);

	auto& vec_neighbours = config_t.get_neighbours();
	const auto& vec_pa_t = config_t.get_particle();
	//size_t total_neighbours = 0;
	//std::cout << vec_neighbours.size() << "\n";

	//size_t index;
	//std::cout << "enter a number<4320: \n";
	//std::cin >> index;
	//size_t i_size = vec_neighbours[index].pvec_neighbours.size();
	//auto& pvec_nbr = vec_neighbours[index].pvec_neighbours;
	//std::cout << "particle(i) id: " << vec_neighbours[index].p_center_pa->id << '\n';
	//std::cout << "neighbour list size: " << i_size << "\n";
	//cout << D2::DELTAij;
	vector<double> vecd_D2;
	vecd_D2.resize(vec_pa_t.size());
	for (size_t i = 0; i < vec_pa_t.size(); i++)
	{
		const Particle& pa = config_t.get_particle()[i];
		D2::D2 d2_pa(&pa, r_cut, &config_t, &config_t_);
		vecd_D2[i] = d2_pa.get_D2();
	}
	//cout << vecd_D2[0].get
	config_t.para_to_dump("D2_test", { "D2" }, { vecd_D2 });
	//D2::d2_pa()
		/*
		std::cout << "particle j id: " << pvec_nbr[0]->id << '\n';
		const Particle& pai = *vec_neighbours[index].p_center_pa;
		const Particle& paj = config_t.get_particle(pvec_nbr[0]->id);
		auto rij = config_t.get_PBC_rij(pai, paj);
		double drij = sqrt(rij.dot(rij));
		cout << "rij:" << '\n' << rij << '\n';
		cout << "drij: " << drij << '\n';
		cout << "i " << pai.id << ": " << pai.rx << ' ' << pai.ry << ' ' << pai.rz << '\n';
		cout << "j " << paj.id << ": " << paj.rx << ' ' << paj.ry << ' ' << paj.rz << '\n';
		auto ri_box = config_t.to_box_coordination(pai);
		auto rj_box = config_t.to_box_coordination(paj);
		cout << "ri box: \n" << ri_box << '\n';
		cout << "rj box: \n" << rj_box << '\n';
		D2::VectorDd rij_box = ri_box - rj_box;
		for (size_t i = 0; i < rij_box.size(); i++)
		{
			rij_box[i] -= round(rij_box[i]);
		}
		cout << "rij box: \n" << rij_box << '\n';
		auto rij_car = config_t.get_m_vector_box_to_cartesian() * rij_box;
		cout << "m_vector_box_to_cartesian * rij_box: \n" << rij_car << '\n';

		cout << config_t.get_m_vector_box_to_cartesian() << '\n';
		*/
	return 0;
}
