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
#include "chrono/core/ChFileutils.h"
#include "SystemGeometry.h"
#include <memory>
namespace chrono {
	namespace irrlicht {
		//namespace irr{
		//	namespace gui{

		class IrrGui : public irr::IEventReceiver {
		public:
			IrrGui(ChIrrApp* myapp, std::vector<std::shared_ptr<Smarticle>> *mySmarticlesVec, std::shared_ptr<SystemGeometry> msys);
			bool OnEvent(const irr::SEvent& event);
			void drawSmarticleAmt(int numLayers);
			void drawAngle();
			void drawSuccessful();
			void drawSuccessful2();
			void drawCamera();
			void addSuccessful(Smarticle &sPtr);
			void resetSuccessfulCount();
			void GenerateVibrateGait(Smarticle& sPtr);
			void SaveToMovie();
			void DeleteImgs();
			void screenshot(int stepsPerFrame);
			int fps = 30;
			int dtPerFrame = 134;
		private:
			std::vector<std::shared_ptr<Smarticle>> *sv;
			ChIrrApp* app;
			std::shared_ptr<SystemGeometry> sys;
			irr::gui::IGUIScrollBar* scrollbar_friction;
			irr::gui::IGUIStaticText* text_Q;
			irr::gui::IGUIStaticText* text_W;
			irr::gui::IGUIStaticText* text_E;
			irr::gui::IGUIStaticText* text_R;
			irr::gui::IGUIStaticText* text_T;
			irr::gui::IGUIStaticText* text_Y;
			irr::gui::IGUIStaticText* text_SmarticleAmt;
			irr::gui::IGUIStaticText* text_Angle;
			irr::gui::IGUIStaticText* text_cameraPos;
			irr::gui::IGUIScrollBar* scrollbar_cohesion;
			irr::gui::IGUIStaticText* text_cohesion;
			irr::gui::IGUIScrollBar* scrollbar_compliance;
			irr::gui::IGUIStaticText* text_compliance;
			irr::gui::IGUIStaticText* text_angle1;
			irr::gui::IGUIStaticText* text_angle2;
			irr::gui::IGUIEditBox* angle1Input;
			irr::gui::IGUIEditBox* angle2Input;
			irr::gui::IGUIEditBox* pgainInput;
			irr::gui::IGUIEditBox* igainInput;
			irr::gui::IGUIEditBox* dgainInput;

			irr::gui::IGUIStaticText* text_successful;
			static int successfulCount;
		};

	};
};
#endif /* GUI_H_ */