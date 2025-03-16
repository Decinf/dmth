#pragma once
#include "dmth.hpp"
#include "SFML/Graphics.hpp"

#define grapth_panel_offset 50
#define char_size 11

namespace dmth
{   
    sf::Vector2f vec2_to_sfvec2(vec2 v){return sf::Vector2f(v.x,v.y);}
    sf::Vector2i vec2_to_sfvec2i(vec2 v){return sf::Vector2i(v.x,v.y);}


    class graph_line_2d
    {
    public:
        sf::Color clr = sf::Color::Red;
        sf::Color point_clr = sf::Color::Red;
        float thicc = 1.0f;
        bool sign_y=false, sign_x=false;
        float point_radius = 0.0f;
        std::vector<vec2> values;
        graph_line_2d(sf::Color clr, float thicc, std::vector<vec2> values){
            this->clr = clr;
            this->thicc = thicc;
            this->values = values;
            this->point_clr = clr;
        }

        graph_line_2d(sf::Color clr, float thicc){
            this->clr = clr;
            this->thicc = thicc;
            this->values = std::vector<vec2>();
            this->point_clr = clr;
        }
        virtual void update_values(float start, float end, float step = 0.1f){};
        
    };

    class graph_line_2d_function : public graph_line_2d
    {
        float (*f)(float);
        float (*fp)(float,float*)=nullptr;
        
    public:
        graph_line_2d_function(sf::Color clr, float thicc, float (*f)(float)) : graph_line_2d(clr,thicc,std::vector<vec2>{}){
            this->f = f;
        }
        graph_line_2d_function(sf::Color clr, float thicc, float (*fp)(float,float*), float*params) : graph_line_2d(clr,thicc,std::vector<vec2>{}){
            this->fp = fp;
            this->params = params;
        }

        float* params = nullptr;

        virtual void update_values(float start, float end, float step = 0.1f){
            values.clear();
            if(params == nullptr)
                for(float x = start; x <= end; x+=step){
                    values.push_back(vec2(x,f(x)));                
                }
            else
                for(float x = start; x <= end; x+=step){
                    values.push_back(vec2(x,fp(x,params)));                
                }
        }
        
    };
    
    class graph_shape_2d
    {
    public:
        sf::Color outline = sf::Color::White;
        sf::Color bg = sf::Color::Transparent;
                
        std::vector<vec2> values;
        graph_shape_2d(sf::Color outline, sf::Color bg, std::vector<vec2> values){
            this->outline = outline;
            this->bg = bg;
            this->values = values;
        }

        graph_shape_2d(sf::Color outline, sf::Color bg){
            this->outline = outline;
            this->bg = bg;
            this->values = std::vector<vec2>();
        } 
        graph_shape_2d(){this->values = std::vector<vec2>();}         
    };



    
    class graph_2d
    {
    private:
        std::string ftostr(float f, uint8_t pop_symbols = 2){
            std::string ret = std::to_string(f);
            for(uint8_t i = 0; i < pop_symbols; ++i)ret.pop_back();
            return ret;
        };

    public:
        std::vector<graph_line_2d*> lines;
        std::vector<graph_shape_2d*> shapes;


        graph_2d(){lines = std::vector<graph_line_2d*>{}; shapes = std::vector<graph_shape_2d*>{};};
        graph_2d(float scale, float grid_scale){this->scale = vec2(scale,scale); this->grid_scale=vec2(grid_scale,grid_scale); lines = std::vector<graph_line_2d*>{}; shapes = std::vector<graph_shape_2d*>{};};
        graph_2d(vec2 scale, vec2 grid_scale){this->scale = scale; this->grid_scale=grid_scale; lines = std::vector<graph_line_2d*>{}; shapes = std::vector<graph_shape_2d*>{};};
        ~graph_2d(){};
        uint wx = 500, wy = 500;
        
        vec2 scale = vec2(1.0f,1.0f); //real graph scaling
        vec2 offset = vec2(0.0f,0.0f); //real graph offsets
        vec2 grid_scale = vec2(1.0f,1.0f);

        sf::Color graph_color = sf::Color::Red;
		sf::Color zero_color = sf::Color(255, 255, 255, 64);
		sf::Color grid_color = sf::Color(255, 255, 255, 48);
		sf::Color bg_color = sf::Color::Black;
        sf::Color fg_color = sf::Color::White;
        sf::Clipboard inp_text = sf::Clipboard();

        vec2 rp(vec2 v){return vec2(v.x/scale.x*wx, -v.y/scale.y*wy);} //real-pixel (coordinates)
        vec2 pr(vec2 v){return vec2(v.x*scale.x/wx, -v.y*scale.y/wy);} //pixel-real (coordinates)

        float rpx(float v){return v/scale.x*wx;} //real-pixel (coordinates)
        float prx(float v){return v*scale.x/wx;} //pixel-real (coordinates)

        float rpy(float v){return -v/scale.y*wy;} //real-pixel (coordinates)
        float pry(float v){return -v*scale.y/wy;} //pixel-real (coordinates)

        sf::Font font;
        bool loadfont(std::string path){return (bool)font.loadFromFile(path);}


        void set_foreground_color(sf::Color color){
            fg_color = color;
            zero_color = fg_color; zero_color.a *= 0.25f;
            grid_color = fg_color; grid_color.a *= 0.2f;
        };

        void set_background_color(sf::Color color){
            bg_color = color;
        };

        int run(std::string wname, std::string screen_save_path="", bool close_instantly=false){
            if(!loadfont("arial.ttf")) return 1;

            sf::Text cord_sign;
            cord_sign.setFillColor(fg_color);
            cord_sign.setCharacterSize(11);
            cord_sign.setFont(font);


            sf::RenderWindow window(sf::VideoMode(wx+grapth_panel_offset, wy+50), wname);
            window.setFramerateLimit(12);

            sf::RectangleShape ox_axis(sf::Vector2f(wx, 3));
            sf::RectangleShape oy_axis(sf::Vector2f(3, wy));

            sf::RectangleShape ox_grid(sf::Vector2f(wx, 1));
            sf::RectangleShape oy_grid(sf::Vector2f(1, wy));

            sf::RectangleShape graph_box(sf::Vector2f(wx, wy));
            graph_box.setFillColor(sf::Color::Transparent);
            graph_box.setOutlineColor(fg_color);
            graph_box.setOutlineThickness(1);

            
            ox_axis.setFillColor(zero_color);
            oy_axis.setFillColor(zero_color);

            ox_grid.setFillColor(grid_color);
            oy_grid.setFillColor(grid_color);

            vec2 dir;
            for(int i = 0; i < lines.size(); ++i){
                lines.at(i)->update_values(offset.x,offset.x+scale.x,scale.x*0.01f);
            }
            while (window.isOpen())
            {
                window.clear(bg_color);
                sf::Event event;
                
                dir=vec2(0,0);
                //std::cout << "afidfngjkfdngk" << std::endl;
                dir.x= (float)sf::Keyboard::isKeyPressed(sf::Keyboard::Right)*grid_scale.x;
                dir.x-= (float)sf::Keyboard::isKeyPressed(sf::Keyboard::Left)*grid_scale.x;
                dir.y= (float)sf::Keyboard::isKeyPressed(sf::Keyboard::Up)*grid_scale.y;
                dir.y-= (float)sf::Keyboard::isKeyPressed(sf::Keyboard::Down)*grid_scale.y;
                /*update values if offset changed*/
                offset = offset + dir;
                if(dir.x != 0.0f || dir.y != 0.0f)
                    for(int i = 0; i < lines.size(); ++i){
                        lines.at(i)->update_values(offset.x,offset.x+scale.x,scale.x*0.01f);
                    }
                vec2 pof = rp(offset); //pixel offset
                
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Q)) {scale = scale * 10.0f; grid_scale = grid_scale * 10.0f;}
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::E)) {scale = scale * 0.1f; grid_scale = grid_scale * 0.1f;}

                
                window.draw(ox_axis);
                window.draw(oy_axis);
                
                std::string numbuff;
                cord_sign.setRotation(90);

                ox_axis.setPosition(0, -pof.y);
                oy_axis.setPosition(-pof.x, 0);

                //create grids in loop
                for (float p = 0.0f; p < scale.x; p+=grid_scale.x)
                {
                    cord_sign.setString(ftostr(p+offset.x));
                    oy_grid.setPosition(rpx(p), 0);
                    cord_sign.setPosition(rpx(p)+8 , wy+char_size);
                    window.draw(oy_grid);
                    window.draw(cord_sign);
                }
                cord_sign.setRotation(0);
                for (float p = -scale.y; p < 0.0f; p+=grid_scale.y)
                {
                    cord_sign.setString(ftostr(p+offset.y));
                    ox_grid.setPosition(0, rpy(p));
                    cord_sign.setPosition(wx+char_size , rpy(p)-8);
                    window.draw(ox_grid);
                    window.draw(cord_sign);
                }

                for(const graph_line_2d * line : lines){
                    sf::CircleShape point;
                    point.setRadius(line->point_radius);
                    point.setFillColor(line->point_clr);
                    if(line->values.size() < 2) continue;
                    for(int i = 0; i < line->values.size()-1;++i){
                        vec2 start = line->values.at(i); start = start - offset;
                        vec2 end = line->values.at(i+1); end = end - offset;
                        if(line->sign_y){
                            sf::Vector2f signpos = vec2_to_sfvec2(rp(start));
                            signpos.y -= 11;
                            cord_sign.setPosition(signpos);
                            cord_sign.setString(ftostr(line->values.at(i).y,3));
                            window.draw(cord_sign);

                            if (i == line->values.size()-2){
                                signpos = vec2_to_sfvec2(rp(end));
                                signpos.y -= 11;
                                cord_sign.setPosition(signpos);
                                cord_sign.setString(ftostr(line->values.at(i+1).y,3));
                                window.draw(cord_sign);
                            }
                        }

                        if(line->point_radius > 0.01f){
                            sf::Vector2f signpos = vec2_to_sfvec2(rp(start));
                            //signpos.y -= 11;
                            point.setPosition(signpos);
                            window.draw(point);

                            if (i == line->values.size()-2){
                                signpos = vec2_to_sfvec2(rp(end));
                                point.setPosition(signpos);
                                window.draw(point);
                            }
                        }
                        
                        sf::Vertex poly_[] =
                        {
                            sf::Vertex(vec2_to_sfvec2(rp(start)),line->clr),
                            sf::Vertex(vec2_to_sfvec2(rp(end)),line->clr)
                        };
                        window.draw(poly_, 2, sf::Lines);
                    }
                }

                for(const graph_shape_2d * shape : shapes){
                    if(shape->values.size() < 3) continue;
                    sf::ConvexShape convex; convex.setPointCount(shape->values.size());
                    convex.setFillColor(shape->bg);
                    convex.setOutlineColor(shape->outline);
                    convex.setOutlineThickness(1);
                    for(int i = 0; i < shape->values.size();++i){
                        vec2 p = shape->values.at(i); p = rp(p - offset);
                        convex.setPoint(i, sf::Vector2f(p.x, p.y));
                    }
                    window.draw(convex);
                }

                /*for (int p = 0; p < poly_count; ++p)
                {
                    float x_s = (-0.5 * ww * scale) + (p * ww * scale)/poly_count;
                    float x_e = (-0.5 * ww * scale) + ((p+1) * ww * scale) / poly_count;
                    float y_s = y0_ + f(x_s-x0_);
                    float y_e = y0_ + f(x_e - x0_);



                    sf::Vertex poly_[] =
                    {
                        sf::Vertex(sf::Vector2f(rxtp(x_s), rytp(y_s))),
                        sf::Vertex(sf::Vector2f(rxtp(x_e), rytp(y_e)))
                    };
                    
                    window.draw(poly_, 2, sf::Lines);
                }*/
                while (window.pollEvent(event))
                {
                    // "close requested" event: we close the window
                    if (event.type == sf::Event::Closed || close_instantly)
                    {
                        if(screen_save_path != ""){
                            sf::Texture texture;
                            texture.create(window.getSize().x, window.getSize().y);
                            texture.update(window);
                            if (texture.copyToImage().saveToFile(screen_save_path))
                            {
                                std::cout << "screenshot saved to " << screen_save_path << std::endl;
                            }
                        }
                        window.clear();
                        window.close();
                        return 0;
                    }
                }
                window.draw(graph_box);
                //window.draw(rs);
                window.display();

            }


        }
        
    };

}