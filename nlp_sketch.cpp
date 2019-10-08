struct NLP {

	struct Objective : VectorFunction {

		inline void compute();

		inline void jacobian();

		VectorFunction* fun;

	};

	Objective obj;

	//--------------------------------------------------------------------------

	struct EqualityConstraints : VectorFunction {

		EqualityConstraints(){
			nEq = 0;
			computeIdx.push_back(0);
		}

		template<class InType, class OutType>
		inline void compute(const MatrixBase<InType> & x, MatrixBase<OutType> const & fx_) const {
			
			MatrixBase<OutType> & fx = fx_.const_cast_derived();

			EqThreads.restart();
			for(int i=0; i<nEq; i++){
				EqThreads.push( [this, x, fx](int thr){ this->funx[i].compute(x, fx.segment(...); ) } );
			}
			EqThreads.stop();

		}

		template<class InType, class JacType>
		inline void jacobian(const MatrixBase<InType> & x, MatrixBase<JacType> const & jx_) const {

			MatrixBase<JacType> jx = jx_.const_cast_derived();

			EqThreads.restart();
			for(int i=0; i<nEq; i++){
				EqThreads.push( [this, x, jx](int thr){ this->funx[i].jacobian(x, jx.block(...)); } );
			}

		}

		inline void addEqConst(VectorFunction* fun){
			this->funx.push_back(fun);
			this->computeIdx.push_back( computeIdx.back() + fun->OR() );
			this->nEq++;
		}

		std::vector<VectorFunction> funx;
		std::vector<int> computeIdx;  // computeIdx[i] = location where funciton 'i' puts data in fx
		int nEq;

		ctpl::thread_pool EqThreads;

	};

	EqualityConstraints eqConsts;

	inline void addEqConst(VectorFunction* fun) {
		this->eqConsts.addEqConst(fun);;
	}

	//--------------------------------------------------------------------------

	struct InequalityConstraints : VectorFunction {

		inline void compute();

		inline void jacobian();

		std::vector<VectorFunciton> funx;

		ctpl::thread_pool InThreads;

	};

};


template<class NLP>
struct NLP_Solver {
	NLP* nlp;
};