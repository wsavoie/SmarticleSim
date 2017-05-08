#include "IrrGui.h"
#include <stdlib.h>  // system, rand, srand, RAND_MAX

using namespace chrono;
using namespace chrono::irrlicht;
using namespace irr;
using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;


IrrGui::IrrGui(ChIrrApp* myapp, std::vector<std::shared_ptr<Smarticle>> *mySmarticlesVec, std::shared_ptr<SystemGeometry>  msys) {

	sv = mySmarticlesVec;
	sys = msys;
	// store pointer applicaiton
	app = myapp;
	// ..add a GUI slider to control friction


	if (bucketType == STRESSSTICK || bucketType == CYLINDER)
	{
		text_Angle = app->GetIGUIEnvironment()->addStaticText(L"PID CONTROL",
			rect<s32>(850, 45, 1050, 60), true);
	}
	else
	{
		text_Angle = app->GetIGUIEnvironment()->addStaticText(L"Angle: 0, Increment: 0",
			rect<s32>(850, 45, 1050, 60), true);
	}
	text_cameraPos = app->GetIGUIEnvironment()->addStaticText(L"Camera Pos and Target",
		rect<s32>(850, 25, 1150, 40),true);
		text_SmarticleAmt = app->GetIGUIEnvironment()->addStaticText(L"Layers: 0, Smarticles: 0",
			rect<s32>(850, 65, 1050, 80), true);

		text_Q = app->GetIGUIEnvironment()->addStaticText(L"Press Q to activate GUI1",
			rect<s32>(850, 85, 1050, 100), true);
		text_W = app->GetIGUIEnvironment()->addStaticText(L"Press W to activate GUI2",
			rect<s32>(850, 105, 1050, 120), true);
		text_E = app->GetIGUIEnvironment()->addStaticText(L"Press E to activate GUI3",
			rect<s32>(850, 125, 1050, 140), true);
		text_R = app->GetIGUIEnvironment()->addStaticText(L"Press R to vibrate around current angle",
			rect<s32>(850, 145, 1050, 160), true);
		text_T = app->GetIGUIEnvironment()->addStaticText(L"Press T to vibrate around specified angles",
			rect<s32>(850, 165, 1050, 180), true);

		angle1Input = app->GetIGUIEnvironment()->addEditBox(L"75",
			rect<s32>(1050, 165, 1100, 180), true);
		angle2Input = app->GetIGUIEnvironment()->addEditBox(L"75",
			rect<s32>(1100, 165, 1150, 180), true);

		text_Y = app->GetIGUIEnvironment()->addStaticText(L"Press Y to move cylinder away",
			rect<s32>(850, 185, 1050, 200), true);
		text_successful = app->GetIGUIEnvironment()->addStaticText(L"Successfully moved smarticles",
			rect<s32>(850, 205, 1050, 220), true);
		resetSuccessfulCount();
		
		wchar_t pm[100]; swprintf(pm,100, L"%g", p_gain);
		wchar_t im[100]; swprintf(im,100, L"%g", i_gain);
		wchar_t dm[100]; swprintf(dm,100, L"%g", d_gain);
		pgainInput = app->GetIGUIEnvironment()->addEditBox(pm,
			rect<s32>(1050, 225, 1100, 240), true);
		igainInput = app->GetIGUIEnvironment()->addEditBox(im,
			rect<s32>(1100, 225, 1150, 240), true);
		dgainInput = app->GetIGUIEnvironment()->addEditBox(dm,
			rect<s32>(1150, 225, 1200, 240), true);
		//std::to_wstring(d_gain).c_str(),
		pgainInput->setID(1000);
		igainInput->setID(1001);
		dgainInput->setID(1002);

		inc = .1;
		//box_ang = 10 * D2R;
		//rampInc = 1.0 / 60.0;
		drum_omega = 1 * 2 * PI;
	}
	int IrrGui::successfulCount = 0;
	bool IrrGui::OnEvent(const SEvent& event) {
		// check if user moved the sliders with mouse..
		if (event.EventType == irr::EET_KEY_INPUT_EVENT && !event.KeyInput.PressedDown) 
		{

			switch (event.KeyInput.Key)
			{
			case irr::KEY_F1:
			{
				sv->at(0)->ChangeActive(!sv->at(0)->active);
				return true;
				break;
			}
			case irr::KEY_F2:
				sv->at(1)->ChangeActive(!sv->at(1)->active);
				return true;
				break;
			case irr::KEY_F3:
				sv->at(2)->ChangeActive(!sv->at(2)->active);
				return true;
				break;
			case irr::KEY_F4:
				sv->at(3)->ChangeActive(!sv->at(3)->active);
				return true;
				break;
			case irr::KEY_F5:
				sv->at(4)->ChangeActive(!sv->at(4)->active);
				return true;
				break;

			case irr::KEY_KEY_Q:
				if (Smarticle::global_GUI_value != MoveType::GUI1)
					Smarticle::global_GUI_value = MoveType::GUI1;
				else
					Smarticle::global_GUI_value = MoveType::GLOBAL;
				return true;
				break;

			case irr::KEY_KEY_W:
				if (Smarticle::global_GUI_value != MoveType::GUI2)
					Smarticle::global_GUI_value = MoveType::GUI2;
				else
					Smarticle::global_GUI_value = MoveType::GLOBAL;
				return true;
				break;
			case irr::KEY_KEY_E:
				if (Smarticle::global_GUI_value != MoveType::GUI3)
					Smarticle::global_GUI_value = MoveType::GUI3;
				else
					Smarticle::global_GUI_value = MoveType::GLOBAL;
				return true;
				break;
			case irr::KEY_KEY_9:
				GetLog() << "Saving frames";
				saveFrame = true;
				
				break;
			case irr::KEY_KEY_8:
			{
				app->GetDevice()->closeDevice();
				break;
			}
			case irr::KEY_KEY_N:
				if (Smarticle::global_GUI_value != MoveType::MIDT)
					Smarticle::global_GUI_value = MoveType::MIDT;
				else
					Smarticle::global_GUI_value = MoveType::GLOBAL;
				return true;
				break;
			case irr::KEY_KEY_R: //vibrate around current theta
				if (Smarticle::global_GUI_value != MoveType::VIB)
				{

					double CurrTheta01;
					double CurrTheta12;
					Smarticle::global_GUI_value = MoveType::VIB;
					std::pair<double, double> angPair;
					for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
					{

						std::shared_ptr<Smarticle> sPtr = sv->at(i);

						CurrTheta01 = sPtr->GetAngle(0);
						CurrTheta12 = sPtr->GetAngle(1);
						sPtr->AssignState(MoveType::VIB);
						sPtr->vib.clear(); //since mv points to address containing current movetype mv clears vib
						sPtr->GenerateVib(CurrTheta01, CurrTheta12);
						//sPtr->addInterpolatedPathToVector(CurrTheta01,CurrTheta12, CurrTheta01 + vibAmp, CurrTheta12 + vibAmp); //curr				->		curr+vib
						//sPtr->addInterpolatedPathToVector(CurrTheta01 + vibAmp, CurrTheta12 + vibAmp, CurrTheta01, CurrTheta12);//curr+vib		->		curr
						//sPtr->addInterpolatedPathToVector(CurrTheta01, CurrTheta12, CurrTheta01 - vibAmp, CurrTheta12 - vibAmp);//curr				->		curr-vib
						//sPtr->addInterpolatedPathToVector(CurrTheta01 - vibAmp, CurrTheta12 - vibAmp, CurrTheta01, CurrTheta12);//curr-vib		->		curr
						//sPtr->vib.clear();
						//sPtr->vib.emplace_back(CurrTheta01, CurrTheta12);
						//sPtr->vib.emplace_back(CurrTheta01 - vibAmp, CurrTheta12 - vibAmp);
						//sPtr->vib.emplace_back(CurrTheta01, CurrTheta12);
						//sPtr->vib.emplace_back(CurrTheta01 + vibAmp, CurrTheta12 + vibAmp);
					}
				}

				else
					Smarticle::global_GUI_value = MoveType::GLOBAL;
				return true;
				break;

			case irr::KEY_KEY_T:
				if (Smarticle::global_GUI_value != MoveType::VIB)
				{
					std::pair<double, double> angPair;
					double ang1;
					double ang2;
					Smarticle::global_GUI_value = MoveType::VIB;
					for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
					{

						std::shared_ptr<Smarticle> sPtr = sv->at(i);
						sPtr->AssignState(MoveType::VIB);
						ang1 = wcstod(angle1Input->getText(), NULL)*D2R;
						ang2 = wcstod(angle2Input->getText(), NULL)*D2R;
						sPtr->vib.clear();

						//in case strange values are written
						if (ang1 > PI || ang1 < -PI)
						{
							Smarticle::global_GUI_value = MoveType::GLOBAL;
							return true;
							break;
						}
						if (ang2 > PI || ang2 < -PI)
						{
							Smarticle::global_GUI_value = MoveType::GLOBAL;
							return true;
							break;
						}
						sPtr->GenerateVib(ang1, ang2);
						//sPtr->addInterpolatedPathToVector(ang1, ang2, ang1 + vibAmp, ang2 + vibAmp);//curr				->		curr+vib
						//sPtr->addInterpolatedPathToVector(ang1 + vibAmp, ang2 + vibAmp, ang1, ang2);//curr+vib		->		curr
						//sPtr->addInterpolatedPathToVector(ang1, ang2, ang1 - vibAmp, ang2 - vibAmp);//curr				->		curr-vib
						//sPtr->addInterpolatedPathToVector(ang1 - vibAmp, ang2 - vibAmp, ang1, ang2);//curr-vib		->		curr

					}
				}
				else
					Smarticle::global_GUI_value = MoveType::GLOBAL;
				return true;
				break;

			case irr::KEY_KEY_Y: //remove container or floor
				if (bucket_exist)
				{
					switch (bucketType)
					{
					case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
						for (size_t i = 0; i < bucket_bod_vec.size(); i++)
						{
							bucket_bod_vec.at(i)->SetBodyFixed(false);

							bucket_bod_vec.at(i)->SetPos(ChVector<>(
								bucket_bod_vec.at(i)->GetPos().x(),
								bucket_bod_vec.at(i)->GetPos().y() + 1,
								bucket_bod_vec.at(i)->GetPos().z()));
						}
						break;
					case HOPPER:
						sys->bucket_bott->SetPos(ChVector<>(1, 0, 0));
						break;
					case FLATHOPPER:
						sys->bucket_bott->SetPos(ChVector<>(1, 4, 0));
						break;
					}

					bucket_exist = false;
				}

				return true;
				break;
			case irr::KEY_KEY_J:
			{
				static bool ran = false;
				switch (bucketType)
				{
				case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
					for (size_t i = 0; i < bucket_bod_vec.size(); i++)
					{
						//bucket_bod_vec.at(i)->SetBodyFixed(false);
						bucket_bod_vec.at(i)->SetBodyFixed(true);

						double x = bucket_bod_vec.at(i)->GetPos().x();
						double y = bucket_bod_vec.at(i)->GetPos().y();
						double z = bucket_bod_vec.at(i)->GetPos().z();
						//GetLog() << "\n" << bucket_bod_vec.at(i).get_ptr()->GetPos() << "\n";

						double theta = atan2(y, x);
						//double theta = atan(x/y);
						double r = sqrt(x*x + y*y);
						bucket_bod_vec.at(i)->SetPos(ChVector<>(x - sys->bucket_rad*.01*cos(theta), y - sys->bucket_rad*.01*sin(theta), z));

					}

					break;
				case HOPPER:
					sys->bucket_bott->SetPos(ChVector<>(100, 0, 0));
					break;
				case FLATHOPPER:
					sys->bucket_bott->SetPos(ChVector<>(100, 0, 0));
					break;
				}

				//return true;
				break;

			}
			case irr::KEY_KEY_K:
			{
				static bool ran = false;
				switch (bucketType)
				{
				case CYLINDER: case STRESSSTICK: case HOOKRAISE: case KNOBCYLINDER:
					for (size_t i = 0; i < bucket_bod_vec.size(); i++)
					{
						//bucket_bod_vec.at(i)->SetBodyFixed(false);
						bucket_bod_vec.at(i)->SetBodyFixed(true);

						double x = bucket_bod_vec.at(i)->GetPos().x();
						double y = bucket_bod_vec.at(i)->GetPos().y();
						double z = bucket_bod_vec.at(i)->GetPos().z();
						//GetLog() << "\n" << bucket_bod_vec.at(i).get_ptr()->GetPos() << "\n";

						double theta = atan2(y, x);
						double r = sqrt(x*x + y*y);
						bucket_bod_vec.at(i)->SetPos(ChVector<>(x + sys->bucket_rad*.01*cos(theta), y + sys->bucket_rad*.01*sin(theta), z));

					}

					break;
				case HOPPER:
					sys->bucket_bott->SetPos(ChVector<>(100, 0, 0));
					break;
				case FLATHOPPER: case BOX:
					sys->bucket_bott->SetPos(ChVector<>(100, 0, 0));
					break;
				}

				//return true;
				break;

			}
			case irr::KEY_KEY_1:
				switch (bucketType)
				{
				case DRUM:
					drum_omega = drum_omega - rampInc*2*PI;
					break;
				case FLATHOPPER: case BOX:
					box_ang = box_ang - rampInc* D2R;
					//rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x() - rampInc * D2R;
					//bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
					//	Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x() - rampInc * D2R
					//	, 0, 0)));
					//bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
					//	Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x() - rampInc * D2R
					//	, 0, 0)));
					break;
				default:
					box_ang = Quat_to_Angle(AngleSet::RXYZ, sys->bucket->GetRot()).x() - rampInc * D2R;
					sys->bucket->SetRot(Angle_to_Quat(AngleSet::RXYZ, ChVector<>(
						Quat_to_Angle(AngleSet::RXYZ, sys->bucket->GetRot()).x() - rampInc * D2R
						, 0, 0)));
					break;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_2:			//increase angle of bucket by rampInc
				switch (bucketType)
				{
				case DRUM:
					drum_omega = drum_omega + rampInc*2*PI;
					break;
				case FLATHOPPER: case BOX:
					box_ang = box_ang + rampInc* D2R;
					//rampAngle = Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x() + rampInc * D2R;
					//bucket->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
					//	Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x() + rampInc * D2R
					//	, 0, 0)));
					//bucket_bott->SetRot(Angle_to_Quat(ANGLESET_RXYZ, ChVector<>(
					//	Quat_to_Angle(ANGLESET_RXYZ, bucket->GetRot()).x() + rampInc * D2R
					//	, 0, 0)));
					break;
				default:
					box_ang = Quat_to_Angle(AngleSet::RXYZ, sys->bucket->GetRot()).x() + rampInc * D2R;
					sys->bucket->SetRot(Angle_to_Quat(AngleSet::RXYZ, ChVector<>(
						Quat_to_Angle(AngleSet::RXYZ, sys->bucket->GetRot()).x() + rampInc * D2R
						, 0, 0)));
					break;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_3:			//decrease rampInc
				rampInc = rampInc - 1;
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_4:			//increase rampInc
				rampInc = rampInc + 1;
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_C:			//decrease p
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					std::shared_ptr<Smarticle> sPtr = sv->at(i);
						p_gain = p_gain - inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_V:			//decrease i
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					std::shared_ptr<Smarticle> sPtr = sv->at(i);
					i_gain = i_gain - inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_B:			//decrease d
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					std::shared_ptr<Smarticle> sPtr = sv->at(i);
					d_gain = d_gain - inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_D:			//increase p
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					std::shared_ptr<Smarticle> sPtr = sv->at(i);
					p_gain = p_gain + inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_F:			//increase i
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					std::shared_ptr<Smarticle> sPtr = sv->at(i);
					i_gain = i_gain + inc;
				}
				drawAngle();
				return true;
				break;
			case irr::KEY_KEY_G:			//increase d
				for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
				{
					std::shared_ptr<Smarticle> sPtr = sv->at(i);
					d_gain = d_gain + inc;
				}
				drawAngle();
				return true;
				break;
			}

		}
		if (event.EventType == EET_GUI_EVENT) 
		{
			s32 id = event.GUIEvent.Caller->getID();
			switch (id)
			{
				case 1000:
					p_gain = wcstod(pgainInput->getText(), NULL);

					break;
				case 1001:
					i_gain = wcstod(igainInput->getText(), NULL);
					break;
				case 1002:
					d_gain = wcstod(dgainInput->getText(), NULL);
					break;
			}
			drawAngle();
		}
		return false;
	}
	void IrrGui::drawSmarticleAmt(int numLayers)//nu
	{
		char message[100]; sprintf(message, "Layers: %d, Smarticles: %d, GUI: %d", numLayers, sv->size(), Smarticle::global_GUI_value);
		this->text_SmarticleAmt->setText(core::stringw(message).c_str());
	}
	void IrrGui::drawAngle()
	{
		if (bucketType == STRESSSTICK || bucketType == CYLINDER)
		{
			if (sv->size() > 0)
			{
				std::shared_ptr<Smarticle> sPtr = sv->at(0);
				char message[100]; sprintf(message, "P:%g, I:%g, D:%g", p_gain, i_gain,d_gain);
				this->text_Angle->setText(core::stringw(message).c_str());
			}
		}
		else if (bucketType == DRUM)
		{
			char message[100]; sprintf(message, "AngVel: %g rpm, Increment: %g", drum_omega/(2*PI) * 60, rampInc * 60);
			this->text_Angle->setText(core::stringw(message).c_str());
		}
		else if (bucketType == BOX)
		{
			char message[100]; sprintf(message, "Angle: %1.3g, Increment: %1.3g", box_ang * R2D, rampInc);
			this->text_Angle->setText(core::stringw(message).c_str());
		}
		else{
			char message[100]; sprintf(message, "Angle: %1.3g, Increment: %1.3g", Quat_to_Angle(AngleSet::RXYZ, sys->bucket->GetRot()).x() * R2D, rampInc);
			this->text_Angle->setText(core::stringw(message).c_str());
		}

	}
	void IrrGui::SaveToMovie()
	{

		char comm[150];
#if defined(_WIN64)
		sprintf(comm, "ffmpeg -framerate %d -i ", fps);
		//strcat(comm, "video_capture/screenshot%05d.png -c:v libxvid -q 0 -vf scale=800:600 video_capture/outVid.avi");
		strcat(comm, "video_capture/screenshot%05d.png -c:v libx264 -crf 17 -vf scale=iw*0.5:ih*0.5 video_capture/outVid.avi" );
#else
		sprintf(comm, "ffmpeg -framerate %d -i ", fps);
		strcat(comm, "video_capture/screenshot%05d.png -c:v libx264 -crf 23 video_capture/outVid.avi" );
#endif
		system(comm);

	}
	void IrrGui::DeleteImgs()
	{
#if defined(_WIN64)
		system("del video_capture\\*.png");
#else
		system("rm -f video_capture/*.png");
#endif
	}
	void IrrGui::screenshot(int stepsPerFrame)
	{
		static int frameNum = 0;
		if (saveFrame) {
			if (frameNum % stepsPerFrame == 0) {
				ChFileutils::MakeDirectory("video_capture");
				#if defined(_WIN64)
				HWND winhandle = reinterpret_cast<HWND>(app->GetVideoDriver()->getExposedVideoData().OpenGLWin32.HWnd);
				//MoveWindow(winhandle, 700, 200, appWidth,appHeight, true);
				SetWindowPos(winhandle, HWND_TOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
				#endif
				irr::video::IImage* image = app->GetVideoDriver()->createScreenShot();
				char filename[100];
				sprintf(filename, "video_capture/screenshot%05d.png", (frameNum + 1) / stepsPerFrame);
				if (image)
					app->GetVideoDriver()->writeImageToFile(image, filename);
				image->drop();
			}
			frameNum++;
		}


		//static int frameNum = 0;
		//double h = app->GetVideoDriver()->getScreenSize().Height;
		//double w = app->GetVideoDriver()->getScreenSize().Width;
		//auto vp = app->GetVideoDriver()->getViewPort();
		////app.GetIGUIEnvironment()->saveGUI("lolol.jpeg",irr::gui::wind);
		//double centx = app->GetVideoDriver()->getViewPort().getCenter().x();
		//double centy = app->GetVideoDriver()->getViewPort().getCenter().y();
		//GetLog() << "Screen Size: " << h << " " << w << "\tScreen: " << centx << " " << centy;
		//if (saveFrame) {
		//	if (frameNum % stepsPerFrame == 0) {
		//		irr::video::IImage* image = app->GetVideoDriver()->createScreenShot();
		//		char filename[100];
		//		sprintf(filename, "video_capture\Myscreenshot%05d.jpeg", (frameNum + 1) / stepsPerFrame);
		//		if (image)
		//			app->GetVideoDriver()->writeImageToFile(image, filename);
		//		image->drop();
		//	}
		//	frameNum++;
		//}
	}
	void IrrGui::drawCamera()
	{
		//application.GetDevice()->getSceneManager()->getRootSceneNode()
		vector3df camPos = app->GetSceneManager()->getActiveCamera()->getAbsolutePosition();
		vector3df camTarget = app->GetSceneManager()->getActiveCamera()->getTarget();
		char message[100]; sprintf(message, "P(%2.3g, %2.3g, %2.3g), T(%2.3g, %2.3g, %2.3g)", camPos.X, camPos.Y, camPos.Z, camTarget.X, camTarget.Y, camTarget.Z);
		this->text_cameraPos->setText(core::stringw(message).c_str());
	}
	void IrrGui::drawSuccessful()
	{
		int count = 0;
		for (size_t i = 0; i < sv->size(); i++) //get each particles current theta
		{
			if (sv->at(i)->successfulMotion) count++;
		}
		char message[100]; sprintf(message, "Successfully Moving: %d/%d", count, sv->size());
		this->text_successful->setText(core::stringw(message).c_str());
	}
	void IrrGui::drawSuccessful2()
	{
		char message[100]; sprintf(message, "Successfully Moving: %d/%d", successfulCount, sv->size());
		this->text_successful->setText(core::stringw(message).c_str());
	}
	void IrrGui::addSuccessful(Smarticle &sPtr)
	{
		if (sPtr.successfulMotion)
			successfulCount++;

	}
	void IrrGui::GenerateVibrateGait(Smarticle& sPtr)
	{
		
	}
	void IrrGui::resetSuccessfulCount()
	{
		successfulCount = 0;
	}