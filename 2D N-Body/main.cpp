#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <unordered_map>

using olc::Key;
using olc::Pixel;

using std::max;
using std::min;
using std::endl;
using std::cout;
using std::vector;
using std::to_string;
using std::unordered_map;

using std::chrono::seconds;
using std::chrono::microseconds;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

#define G 1
#define gravityRange 8
#define controlFriction 0.1

class clusterBall
{
public:
	double x;
	double xv;
	int size;

	clusterBall()
	{
		x = 0;
		xv = 0;
		size = 0;
	}
};

class ball
{
public:
	double x;
	double xv;
	Pixel color;

	ball(double X, double XV, Pixel Color)
	{
		x = X;
		xv = XV;
		color = Color;
	}
};

class Example : public olc::PixelGameEngine
{
public:
	Example() { sAppName = "2D-NBody"; }

	double x;
	double xv;
	double zoom;

	int halfScreen;
	unsigned int m_z;
	unsigned int m_w;

	vector<ball> balls;

	unordered_map<uint32_t, vector<int>> collisionSpace;
	unordered_map<uint32_t, clusterBall> gravitySpace[gravityRange];

	unsigned int intRand()
	{
		m_z = 36969 * (m_z & 65535) + (m_z >> 16);
		m_w = 18000 * (m_w & 65535) + (m_w >> 16);

		return (m_z << 16) + m_w;
	}

	double doubleRand() { return (intRand() + 1.0) * 2.328306435454494e-10; }

	Pixel mapToRainbow(double d)
	{
		d = d - 6.0 * int(d / 6);

		double r = (d > 3) ? max(0.0, min(1.0, d - 4)) : max(0.0, min(1.0, 2 - d));
		double g = (d > 2) ? max(0.0, min(1.0, 4 - d)) : max(0.0, min(1.0, d));
		double b = (d > 4) ? max(0.0, min(1.0, 6 - d)) : max(0.0, min(1.0, d - 2));

		return Pixel(r * 0xff, g * 0xff, b * 0xff);
	}

	void control(double fElapsedTime)
	{
		if (GetKey(Key::Q).bHeld) { zoom /= pow(2, fElapsedTime); }
		if (GetKey(Key::E).bHeld) { zoom *= pow(2, fElapsedTime); }

		if (GetKey(Key::A).bHeld || GetKey(Key::LEFT).bHeld) { xv -= 1000 / zoom * fElapsedTime; }
		if (GetKey(Key::D).bHeld || GetKey(Key::RIGHT).bHeld) { xv += 1000 / zoom * fElapsedTime; }

		xv *= pow(controlFriction, fElapsedTime);
		x += xv * fElapsedTime;
	}

	void ballPullBall(clusterBall* ball1, clusterBall* ball2)
	{
		double dpos = ball2->x - ball1->x;
		double dis = abs(dpos);
		dpos /= dis;
		dpos *= ((long long)ball1->size * ball2->size) * G / (dis * dis);

		ball1->xv += dpos;
		ball2->xv -= dpos;
	}

	void gravity()
	{
		unordered_map<uint32_t, clusterBall>::iterator find1;

		for (int i = 0; i < gravityRange; i++)
		{
			for (unordered_map<uint32_t, clusterBall>::iterator j = gravitySpace[i].begin(); j != gravitySpace[i].end(); j++)
			{
				find1 = gravitySpace[i].find(j->first + 2);

				if (find1 != gravitySpace[i].end())
					ballPullBall(&gravitySpace[i][j->first], &gravitySpace[i][uint32_t(j->first) + uint32_t(2)]);

				if (!(j->first & 1))
				{
					find1 = gravitySpace[i].find(j->first + 3);

					if (find1 != gravitySpace[i].end())
						ballPullBall(&gravitySpace[i][j->first], &gravitySpace[i][uint32_t(j->first) + uint32_t(3)]);
				}

				if (i == 6)
					cout << j->first << ", " << j->first + 2 << " and " << j->first + 3 << endl;
			}
		}

		/*for (unordered_map<uint32_t, clusterBall>::iterator j = gravitySpace[0].begin(); j != gravitySpace[0].end(); j++)
		{
			find1 = gravitySpace[0].find(j->first + 1);

			if (find1 != gravitySpace[0].end())
				ballPullBall(&gravitySpace[0][j->first], &gravitySpace[0][j->first + 1]);
		}*/
	}

	void ballHitBall(ball* ball1, ball* ball2)
	{
		double dpos = ball2->x - ball1->x;
		double dis = abs(dpos);

		if (dis < 1)
		{
			dpos /= dis;
			dis = dpos * (ball2->xv - ball1->xv);

			if (dis < 0)
			{
				dpos *= dis;
				ball1->xv += dpos;
				ball2->xv -= dpos;
			}
		}
	}

	void collision()
	{
		unordered_map<uint32_t, vector<int>>::iterator find1;

		for (unordered_map<uint32_t, vector<int>>::iterator i = collisionSpace.begin(); i != collisionSpace.end(); i++)
		{
			for (int j = 0; j < i->second.size(); j++)
			{
				for (int k = j + 1; k < i->second.size(); k++)
					ballHitBall(&balls[i->second[j]], &balls[i->second[k]]);

				find1 = collisionSpace.find(uint32_t(i->first) + uint32_t(1));

				if (find1 != collisionSpace.end()) {
					for (int k = 0; k < find1->second.size(); k++)
						ballHitBall(&balls[i->second[j]], &balls[find1->second[k]]);
				}
			}
		}
	}

	void drawAndSetSpace(double fElapsedTime)
	{
		Clear(Pixel(0, 0, 0));
		collisionSpace.clear();

		for (int i = gravityRange - 1; i > 0; i--)
		{
			for (unordered_map<uint32_t, clusterBall>::iterator j = gravitySpace[i].begin(); j != gravitySpace[i].end(); j++)
			{
				double mx = (j->second.x - x) * zoom + halfScreen;
				double my = 450 - i * 3 * zoom;
				FillCircle(mx, my, zoom * sqrt(j->second.size) * 0.5);
				DrawLine(mx, my, mx + j->second.xv * 1000, my, Pixel(100, 100, 100));
			}
		}

		for (int i = gravityRange - 1; i > 0; i--)
		{
			for (unordered_map<uint32_t, clusterBall>::iterator j = gravitySpace[i].begin(); j != gravitySpace[i].end(); j++)
			{
				uint32_t bPos = j->first << 1;

				if (gravitySpace[i - 1].find(bPos) != gravitySpace[i - 1].end())
					gravitySpace[i - 1][bPos].xv += j->second.xv;

				bPos ^= 1;

				if (gravitySpace[i - 1].find(bPos) != gravitySpace[i - 1].end())
					gravitySpace[i - 1][bPos].xv += j->second.xv;
			}

			gravitySpace[i].clear();
		}

		for (int i = 0; i < balls.size(); i++)
		{
			balls[i].xv += gravitySpace[0][int(balls[i].x) - (balls[i].x < 0)].xv * fElapsedTime;
			balls[i].x += balls[i].xv * fElapsedTime;
		}

		gravitySpace[0].clear();

		for (int i = 0; i < balls.size(); i++)
		{
			uint32_t floorBPos = int(balls[i].x);
			collisionSpace[floorBPos].push_back(i);
			gravitySpace[0][floorBPos].x += balls[i].x;
			gravitySpace[0][floorBPos].size++;

			FillCircle((balls[i].x - x) * zoom + halfScreen, 450, zoom * 0.5, balls[i].color);
		}

		for (int i = 0; i < gravityRange - 1; i++)
		{
			for (unordered_map<uint32_t, clusterBall>::iterator j = gravitySpace[i].begin(); j != gravitySpace[i].end(); j++)
			{
				j->second.x /= j->second.size;
				uint32_t floorBPos = j->first >> 1;
				gravitySpace[i + 1][floorBPos].x += j->second.x * j->second.size;
				gravitySpace[i + 1][floorBPos].size += j->second.size;
			}
		}

		for (unordered_map<uint32_t, clusterBall>::iterator i = gravitySpace[gravityRange - 1].begin(); i != gravitySpace[gravityRange - 1].end(); i++)
			i->second.x /= i->second.size;
	}

	bool OnUserCreate() override
	{
		x = 0;
		xv = 0;
		zoom = 16;
		halfScreen = (double)ScreenWidth() / 2;
		m_z = (unsigned int)duration_cast<seconds>(high_resolution_clock::now().time_since_epoch()).count();
		m_w = (unsigned int)duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();

		for (int i = 0; i < 100; i++)
		{
			double bPos = doubleRand() * 200 - 100;
			double bPosv = doubleRand() * 0;// 6 - 3;
			Pixel bColor = mapToRainbow(doubleRand() * 6);
			balls.push_back(ball(bPos, bPosv, bColor));
		}
		////balls.push_back(ball(0, 1, mapToRainbow(doubleRand() * 6)));

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		//fElapsedTime *= 10;
		control(fElapsedTime);
		gravity();
		collision();
		drawAndSetSpace(fElapsedTime);

		return true;
	}
};

int main()
{
	/*double nu = -1;
	int num = 0xffffffff;

	for (int i = 0; i < 32; i++)
		cout << (num >> (31 - i) & 1);

	cout << endl << num << endl;

	num = num + 2;

	for (int i = 0; i < 32; i++)
		cout << (num >> (31 - i) & 1);

	cout << endl << num << endl;

	num -= 5;

	for (int i = 0; i < 32; i++)
		cout << (num >> (31 - i) & 1);

	cout << endl << num << endl;

	num /= 2;

	for (int i = 0; i < 32; i++)
		cout << (num >> (31 - i) & 1);

	cout << endl << num << endl;

	return 0;*/

	Example demo;

	if (demo.Construct(1000, 500, 1, 1))
		demo.Start();

	return 0;
}