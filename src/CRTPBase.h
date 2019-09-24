#pragma once
#include "../pch.h"

namespace ASRLib {

	template<class Derived>
	struct CRTPBase {
		Derived& derived() { return static_cast<Derived &>(*this); }
		const Derived& derived() const { return static_cast<const Derived &>(*this); }
		std::string name() const { return std::string(typeid(Derived).name()); }
		std::any asany() const { return std::any(this->derived()); };

		template<class T>
		T erased() const { return T(this->derived()); }

		template<class T>
		T cast() const { return T(this->derived()); }

		template<class T>
		void copyinto(T & obj) const {obj = T(this->derived());}

		template<class T>
		void hardcopyinto(T & obj) const { obj = T(this->derived()); }

		template<class T>
		void copyfrom(T obj) {
			if (this->name() == obj.name()) {
				this->derived() = std::any_cast<Derived>(obj.asany());
			}
			else {
				throw std::invalid_argument("Erased object does not match Derived Type.");
			}
		}
		template<class T>
		Derived unerase(T obj) {
			if (this->name() == obj.name()) {}
			else {
				throw std::invalid_argument("Erased object does not match Derived Type.");
			}
			return std::any_cast<Derived>(obj.asany());
		}
	};
}