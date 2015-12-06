#ifndef IRRGUI_H_
#define IRRGUI_H_

#include <stdlib.h>
#include "common.h"
#include "Smarticle.h"
//#include <irrlicht.h>
//#include "unit_IRRLICHT/ChIrrApp.h"
#include "chrono_irrlicht/ChBodySceneNode.h"  //changed path from unit to chrono to reflect changes in updated chrono
#include "chrono_irrlicht/ChBodySceneNodeTools.h"
//#include "unit_IRRLICHT/ChIrrTools.h"
#include "chrono_irrlicht/ChIrrWizard.h"
#include "core/ChRealtimeStep.h"
//#include <irrlicht.h>
#include "assets/ChTexture.h"
#include <memory>
namespace chrono{
	//namespace irr{
	//	namespace gui{

			class IrrGui : public irr::IEventReceiver {
			public:
				IrrGui(irr::ChIrrApp* myapp, std::vector<Smarticle*> *mySmarticlesVec);
				bool OnEvent(const irr::SEvent& event);
				void drawSmarticleAmt(int numLayers);
				void drawAngle();
				void drawSuccessful();
				void drawSuccessful2();
				void addSuccessful(Smarticle &sPtr);
				void resetSuccessfulCount();


			private:
				std::vector<Smarticle*> *sv;
				irr::ChIrrApp* app;
				irr::gui::IGUIScrollBar* scrollbar_friction;
				irr::gui::IGUIStaticText* text_Q;
				irr::gui::IGUIStaticText* text_W;
				irr::gui::IGUIStaticText* text_E;
				irr::gui::IGUIStaticText* text_R;
				irr::gui::IGUIStaticText* text_T;
				irr::gui::IGUIStaticText* text_Y;
				irr::gui::IGUIStaticText* text_SmarticleAmt;
				irr::gui::IGUIStaticText* text_Angle;
				irr::gui::IGUIScrollBar* scrollbar_cohesion;
				irr::gui::IGUIStaticText* text_cohesion;
				irr::gui::IGUIScrollBar* scrollbar_compliance;
				irr::gui::IGUIStaticText* text_compliance;
				irr::gui::IGUIStaticText* text_angle1;
				irr::gui::IGUIStaticText* text_angle2;
				irr::gui::IGUIEditBox* angle1Input;
				irr::gui::IGUIEditBox* angle2Input;
				irr::gui::IGUIStaticText* text_successful;
				static int successfulCount;
			};

		};
#endif /* GUI_H_ */