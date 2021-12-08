#pragma once
#include "configuration.h"
#include "particle.h"
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <initializer_list>
#include <eigen3/Eigen/Dense>

#define NO_STORE_NEIGHBOUR_RIJ
#define SAFE_GET_NEIGHBOURS
namespace D2
{
	static constexpr size_t DIMENSION = 3; // system dimension, in my case is always 3, other cases NOT tested
	static constexpr size_t NEIGHBOUR_SIZE_REDUNDANCY = 5;//
	using VectorDd = Eigen::Matrix<double, DIMENSION, 1>;// define DIMENSION-d col vector(for position, distance and other concepts in math or physics, NOT a container like std::vector)
	using MatrixDd = Eigen::Matrix<double, DIMENSION, DIMENSION>;// define DIMENSION * DIMENSION shape matrix
	static MatrixDd DELTAij = MatrixDd::Identity();

	struct Neighbour
	{
		const Particle* p_neighbour;
#ifndef NO_STORE_NEIGHBOUR_RIJ //store rij information if necessary
		VectorDd rij;
		double drij;
#endif // !NO_STORE_NEIGHBOUR_RIJ
	};

	struct Neighbours // contains a pointer to center particle and a std::vector container of pointers which points to neighbour-particles of center particle
	{
		const Particle* p_center_pa = nullptr;
		std::vector<Neighbour> vec_neighbours;
		~Neighbours()
		{
			p_center_pa = nullptr;
			vec_neighbours.clear();
		}
	};

	class Configuration_neighbours :public Configuration::Configuration
		/*
		* derived frome base class Configuration::Configuration
		* configuration with neighbour information of all particles
		* nieghbour information stored in std::vector<Neighbours> vec_neighbours
		* these for matrixs are used transfer coordinations from cartesian to box-base, vice versa
		MatrixDd m_base_cartesian_to_box;
		MatrixDd m_base_box_to_cartesian;
		MatrixDd& m_vector_cartesian_to_box = m_base_box_to_cartesian;
		MatrixDd& m_vector_box_to_cartesian = m_base_cartesian_to_box;
		* these transformations are used to compute distance in tilt box under periodic boundary condition, or say PBC
		*/
	{
	public:
		template<typename T>
		T square(T x) { return x * x; }; // define x^2
		void _update_neighbours(); // update neighbours
		Configuration_neighbours(std::string config_file, double _r_cut, BoxType _boxtype = BoxType::orthogonal, PairStyle _pairstyle = PairStyle::single)
			:Configuration(config_file, _boxtype, _pairstyle)
			/*
			* rcut is the cutoff distance of pair-interaction
			* construction
			* update neighbour information
			* setting transfer matrixs
			*/
		{
			rcut = _r_cut;
			m_base_cartesian_to_box = MatrixDd::Zero();

			m_base_cartesian_to_box(0, 0) = get_lx();
			m_base_cartesian_to_box(1, 1) = get_ly();
			m_base_cartesian_to_box(2, 2) = get_lz();
			m_base_cartesian_to_box(0, 1) = get_xy();
			m_base_cartesian_to_box(0, 2) = get_xz();
			m_base_cartesian_to_box(1, 2) = get_yz();
			m_base_box_to_cartesian = m_base_cartesian_to_box.inverse();
			_update_neighbours();
		};

		void sort_neighbours_as_center_pid();
		const std::vector<Neighbours>& get_neighbours() const { return vec_neighbours; }// return corresponding neighbours of ALL particles
		const std::vector<const Neighbours*>& get_sorted_neighbours_pvec() const { return pvec_neighbours_sorted; };
		const Neighbours& get_neighbours(size_t pid) const; // return neighbours of given particle
		const Neighbours& get_neighbours(const Particle& pa) const // return neighbours of given particle
		{
			return get_neighbours(pa.id);
		}

		VectorDd to_box_coordination(const Particle& pa) const; // transfer coordinations of given particle from cartesian to box-base, return a VectorDd
		VectorDd get_PBC_rij(const Particle& pai, const Particle& paj) const; // computing distance between two particles, under PBC

		double get_rcut() { return rcut; } // return rcut

		const MatrixDd& get_m_base_cartesian_to_box() const { return m_base_cartesian_to_box; }
		const MatrixDd& get_m_base_box_to_cartesian() const { return m_base_box_to_cartesian; }
		const MatrixDd& get_m_vector_cartesian_to_box() const { return m_vector_cartesian_to_box; }
		const MatrixDd& get_m_vector_box_to_cartesian() const { return m_vector_box_to_cartesian; }
	private:
		double rcut = 0;
		bool flag_neighbours_update = false;
		bool flag_neighbours_sorted = false;
		std::vector<Neighbours> vec_neighbours;
		std::vector<const Neighbours*> pvec_neighbours_sorted;
		MatrixDd m_base_cartesian_to_box; // matrix, transfer base vectors from cartesian to box-base
		MatrixDd m_base_box_to_cartesian; // matrix, inverse of m_base_cartesian_to_box
		MatrixDd& m_vector_cartesian_to_box = m_base_box_to_cartesian;
		MatrixDd& m_vector_box_to_cartesian = m_base_cartesian_to_box;
	};


	class D2 // compute D2 of given particle
	{
	public:
		D2() {};
		D2(const Particle* _pa, double _rcut, const Configuration_neighbours* _config_t, const Configuration_neighbours* _config_t_)
			/*
			* compute D^2_min, see:
			* Falk, M. L., & Langer, J. S. (1998). Dynamics of viscoplastic deformation in amorphous solids. Physical Review E, 57(6), 7192-7205. doi:DOI 10.1103/PhysRevE.57.7192
			* Cubuk, E. D., Ivancic, R. J. S., Schoenholz, S. S., Strickland, D. J., Basu, A., Davidson, Z. S., . . . Liu, A. J. (2017). Structure-property relationships from universal signatures of plasticity in disordered solids. Science, 358(6366), 1033-1037. doi:doi:10.1126/science.aai8830
			* this quantity need particle coordination & neighbours information in two time moments: t, t-¦¤t(note as t_ in this code)
			* all the information need could be gained from configurtaions of two moments
			* BE CAREFULL: though only const pointers are stored, pay attention to life-cycle of variables
			*/
		{
			p_pa = _pa;
			rcut = _rcut;
			config_t = _config_t;
			config_t_ = _config_t_;
		};
		~D2()
		{
			p_pa = nullptr;
			config_t = nullptr;
			config_t_ = nullptr;
		}
		void compute_Xij();
		/* compute Xij, and set flag_computed_Xij to be true. allowing get Xij using get_Yij()
		* results will be stored in this->Xij, which is a MitrixDd(==Eigen::Matrix<Dimension,Dimension,double> type matrix data
		* note: when compute epsilonij or compute D2, Xij will ALWAYS be re-calculated ignoring flag_computed_Xij
		*/
		void compute_Yij();
		/* compute Yij, and set flag_computed_Yij to be true. allowing get Yij using get_Yij()
		* results will be stored in this->Yij, which is a MitrixDd(==Eigen::Matrix<Dimension,Dimension,double> type matrix data
		* note: when compute epsilonij or compute D2, Yij will ALWAYS be re-calculated ignoring flag_computed_Yij
		*/
		void compute_epsilonij();
		/* note: when compute epsilonij, Xij& Yij will ALWAYS be re-calculated, ignoring flag_computed_Xij or flag_computed_Yij;
		* note: when compute D2, epsilonij will ALWAYS be re-calculated, ignoring flag_coputed_epsilonij
		*       which means Xij, Yij and epsilonij will be re=calculated when compute D2
		* run this function will set flag_computed_epsilonij to be true, allowing get epsilonij using get_epsilonij()
		* results will be stored in this->epsilonij, which is a MitrixDd(==Eigen::Matrix<Dimension,Dimension,double> type matrix data
		*/
		void compute_D2();
		/* compute D2, and set flag_computed_D2 to be true, allowing get D2 using get_D2()
		* results will be stored in this->d2, which is a double type scalar data
		*/
		MatrixDd get_Xij()
		{
			if (!flag_computed_Xij) compute_Xij();
			return Xij;
		}
		MatrixDd get_Yij()
		{
			if (!flag_computed_Yij) compute_Yij();
			return Yij;

		}
		MatrixDd get_epsilonij()
		{
			if (!flag_computed_epsilonij) compute_epsilonij(); // note: when compute epsilonij, Xij & Yij will ALWAYS be re-calcualated, ignoring flag_computed_Xij or flag_computed_Yij;
			return epsilonij;
		}
		double get_D2()
		{
			if (!flag_computed_d2) compute_D2();
			return d2;
		}

		void set_pa(const Particle* _pa) { p_pa = _pa; };
	private:
		const Particle* p_pa = nullptr;
		std::vector<Particle> neighbour;
		double d2 = 0.;
		double rcut = 0.;
		const Configuration_neighbours* config_t = nullptr;
		const Configuration_neighbours* config_t_ = nullptr;
		MatrixDd Xij;
		MatrixDd Yij;
		MatrixDd epsilonij;
		bool flag_computed_Xij = false;
		bool flag_computed_Yij = false;
		bool flag_computed_epsilonij = false;
		bool flag_computed_d2 = false;
	};

}



