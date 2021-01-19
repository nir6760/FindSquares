#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include <iostream>
#include <vector>
#include <string.h>
#include <map>
#include <algorithm>
#include <iomanip>
using namespace std;
using namespace cv;


/*
 * section 1: get an image with few squares
 *
 * All the squares are inside the image, each edge is 1 width size pixel.
 * There is no squares overlapping (and it cant be a square inside other square), and each square is at least distance 1 pixel away from a neighbor square.
 * Each square has 1 greyscale value between 1-255, background is 0.
 * We can assume input is valid.
 * 
 * 1. find and return a list of the squares in the image so, each list
 * contains all the properties of the squares, meaning top left coordinate
 * edge width size and the gray-scale value of the square.
 * 
 * 2. do the same as 1 but the solution should be recursive.
 * 
 * 3. Get an image with 1 diagonal square (not horizontal to axis), find this square (this image must be smaller than section 1 image)
 * and find the best overlapping fit square from section 1. Meaning both of the squares one above the other
 * have the biggest overlapped area.
 * Display this overlapping
 *
 * 
 * Basic idea:
 * travarse the image and find all squares,
 * while avoiding traversing areas where the gurntee is there are no squares over there.
 * return the squares properties.
 * Do the same idea recursivly.
 * find the diagonal image in the second image'and fins his best ovelapping square,
 * meaning when you put one over the other, you have the biggest overlapping area.
 * (It is not allowed to rotate the squares)
 */






// class for representing the square.
//for horizontal square we have the row and col coordinates of the top left corner
//for diagonal square we have the row and col coordinates of the top left corner
//  edge size and gray_scale value of square
class Square {
    /*horizontal     diagonal
     * 2 2 2          0 2 0
     * 2 0 2          2 0 2
     * 2 2 2          0 2 0
     * default corner = top left corner (horizontal)
     * default corner = top corner (diagonal)
     * index starts from 0
     */
public:
    uint8_t row_corner_point;
    uint8_t col_corner_point;
    uint8_t edge_size;
    uint8_t gray_scale_color;
    Square(uint8_t row_val, uint8_t col_val, uint8_t edge_size_val, uint8_t gray_scale_val) :
        row_corner_point(row_val), col_corner_point(col_val), edge_size(edge_size_val), gray_scale_color(gray_scale_val) {}
    friend ostream& operator<<(ostream& os, const Square& square);

};
//printing square - overload operator
ostream& operator<<(ostream& os, const Square& square)
{
    os << "top_left_corner_coordinates(" << +square.row_corner_point << ", " << +square.col_corner_point << ')'
        << ", edge_size:" << +square.edge_size << ", pixel_gray_scale:" << +square.gray_scale_color;
    os << endl;//the '+' is for printing the numeric value of uint_8
    return os;
}
#define BACKROUNDCOLOR 0
//declarations of the main functions for better understanding
vector<Square> findAllSquares(vector<vector<uint8_t>>& img); // section 1
vector<Square> findAllSquaresRecursive(vector<vector<uint8_t>>& img); // section 2
Square* findBestOverlap(Square& diagonalSquare, vector<Square>& squares_vec); // section 3
void visualizeOverLapping(Square& diagonalSquare, Square& horizontalSquare);  //section 3


// struct compare helper for the map we are using for visited squares. if the data is sorted, 
struct classCompForMap {
    bool operator() (const pair<uint8_t, uint8_t>& lhs, const pair<uint8_t, uint8_t>& rhs) const
    {
        return lhs.second < rhs.second;
    }
};

// helper function for recursive traverse image declaration
int traverseImage(vector<vector<uint8_t>>& img, int row_iter,
    int col_iter, map<pair<uint8_t, uint8_t>, uint8_t, classCompForMap>& visited,
    vector<Square>& ans);

// check if the coordinate given by row and col is top right corner of square
bool checkUp(vector<vector<uint8_t>>& img, uint8_t row, uint8_t col) {
    if (row - 1 >= 0) {// in size
        if (img[row - 1][col] != BACKROUNDCOLOR) // it's not  new square
            return false;
        return true;
    }
    return true;
}
// find the edge size when given the coordinate of top right corner of square
uint8_t findEdgeSize(vector<vector<uint8_t>>& img, uint8_t row, uint8_t col) {
    uint8_t edge_size_cnt = 0;
    while (col < img[0].size() && img[row][col] != BACKROUNDCOLOR) {
        ++col;
        ++edge_size_cnt;
    }
    return edge_size_cnt;
}

//return the size of the step we can avoid when meeting a known right edge of square,
// delete this square from visited map if we won't meet him again in our search - map which help us avoiding traverse areas we know
//there will be no squares over there.
uint8_t sizeOfStepAndDeleteVisitIfPossible(int row, int col, map<pair<uint8_t, uint8_t>, uint8_t, classCompForMap>& visited) {
    for (auto map_iter = visited.begin();map_iter != visited.end();map_iter++) {
        if (col == (*map_iter).first.second &&
            (*map_iter).first.first < row && row < (*map_iter).first.first + (*map_iter).second) {// enter when found the properties of the square we met again
            uint8_t step = (*map_iter).second;
            if (row == (*map_iter).first.first + (*map_iter).second) //this is the last time we need this square
                visited.erase(map_iter);
            return step;
        }

    }
    return 0; // should not ger here, because it's enter only on visited square

}


//find all the squares in the image and return list of them - section 1
vector<Square> findAllSquares(vector<vector<uint8_t>>& img) {
    vector<Square> ans;
    map<pair<uint8_t, uint8_t>, uint8_t, classCompForMap> visited;
    int row_iter = 0, col_iter = 0;
    while (row_iter < img.size()) {
        col_iter = 0;
        while (col_iter < img[0].size()) {
            if (img[row_iter][col_iter] != BACKROUNDCOLOR) { // this is part of left edge of square
                if (checkUp(img, row_iter,
                    col_iter)) { // check if it's a new square
                    int curr_edge = findEdgeSize(img, row_iter, col_iter);
                    ans.emplace_back(Square(row_iter, col_iter, curr_edge,
                        img[row_iter][col_iter]));
                    visited[make_pair(row_iter,
                        col_iter)] = curr_edge;
                    col_iter = col_iter + curr_edge + 1; // step and go out of the square, we can avoid this space
                }
                else {
                    int step = sizeOfStepAndDeleteVisitIfPossible(row_iter, col_iter, visited);
                    col_iter = col_iter + step + 1; // step and go out of the square, we can avoid this space
                }
            }
            else
                ++col_iter;
        }
        ++row_iter;
    }
    return ans;
}

//find all the squares in the image and return list of them recursively - section 2
//The recursion is by DFS when we avoid squares we found at the past.
// This function use traverseImage for traversing recursively
// The recursive function is not correct way for big images (for example 256*256) because the recursion stack is becoming too big
// for program to handle
// So the program will throw an exeption
//This soloution only work on small scale images(e.g 50*50 or below)
vector<Square> findAllSquaresRecursive(vector<vector<uint8_t>>& img) {
    vector<Square> ans;
    map<pair<uint8_t, uint8_t>, uint8_t, classCompForMap> visited;
    traverseImage(img, 0, 0, visited, ans);
    return ans;
}

// function for traversing the image recursively
int traverseImage(vector<vector<uint8_t>>& img, int row_iter,
    int col_iter, map<pair<uint8_t, uint8_t>, uint8_t, classCompForMap>& visited,
    vector<Square>& ans) {
    // If the entire column is traversed
    if (col_iter >= img[0].size())
        return 0;
    // If the entire row is traversed
    if (row_iter >= img.size())
        return 1;
    // Print the value of the current
    // cell of the matrix
    if (img[row_iter][col_iter] != BACKROUNDCOLOR) { // this is part of left edge of square
        if (checkUp(img, row_iter,
            col_iter)) { // check if it's a new square
            int curr_edge = findEdgeSize(img, row_iter, col_iter);
            ans.emplace_back(Square(row_iter, col_iter, curr_edge,
                img[row_iter][col_iter]));
            visited[make_pair(row_iter,
                col_iter)] = curr_edge;

            return traverseImage(img, row_iter,
                col_iter + curr_edge + 1,
                visited,
                ans);// step and go out of the square, we can avoid this space
        }
        else {
            int step = sizeOfStepAndDeleteVisitIfPossible(row_iter, col_iter,
                visited);
            return traverseImage(img, row_iter,
                col_iter + step + 1,
                visited,
                ans);// step and go out of the square, we can avoid this space
        }
    }
    // Recursive call to traverse the matrix
    // in the Horizontal direction
    if (traverseImage(img, row_iter,
        col_iter + 1, visited, ans) == 1)
        return 1;
    // Recursive call for changing the
    // Row of the matrix
    return traverseImage(img, row_iter + 1, 0, visited, ans);
}

// check if the coordinate are inside the image and the point is still part of the right top edge of diagonal square
bool insideImg(vector<vector<uint8_t>>& img, int row, int col) {
    return row < img.size() && col < img[0].size()
        && img[row][col] != BACKROUNDCOLOR; // so its valid
}
//exeption for no squares
struct noSquareException : public exception {
    const char* what() const throw () {
        return "There is no squares in imag";
    }
};
// find the diagonal square in the image
//if there is no square, throw an exeption (input is not valid)
Square* findDiagonalSquare(vector<vector<uint8_t>>& img) {
    // I will first meet his top right corner
    for (int row_iter = 0; row_iter < img.size(); ++row_iter) {
        for (int col_iter = 0; col_iter < img[0].size(); ++col_iter) {

            if (img[row_iter][col_iter] != BACKROUNDCOLOR) { // this is the top corner of the diagonal square
                int cnt_size = 0;
                int i = row_iter, j = col_iter;
                // go along the right edge to its end
                while (insideImg(img, i, j)) {
                    ++i;
                    ++j;
                    ++cnt_size;
                }
                return new Square(row_iter, col_iter, cnt_size, img[row_iter][col_iter]);
            }
        }

    }
    throw noSquareException(); // should not get here if the input is valid
}

// struct helper for comparing areas of diagonal square and horizontal square
struct compareAreas {
    int diagonalArea;
    explicit compareAreas(int diagonalArea_val) :diagonalArea(diagonalArea_val) {}
    bool operator() (Square& s1, Square& s2) {
        int distFrom_s1 = abs(s1.edge_size * s1.edge_size - diagonalArea);
        int distFrom_s2 = abs(s2.edge_size * s2.edge_size - diagonalArea);
        return (distFrom_s1 < distFrom_s2);
    }
};


// find the best overlapping square from the vector list for the diagonal square
// the best overlapping one is the one with the similar area (number of pixels inside the square include edges) to the diagonal square
// if there is no best overlapping square (meaning there are no horizontal squares) than throw exception
Square* findBestOverlap(Square& diagonalSquare, vector<Square>& squares_vec) {
    if (squares_vec.empty()) // there are no squares in vec
        throw noSquareException();
    // way for finding area of horizontal square
    int  diagonalArea = diagonalSquare.edge_size * diagonalSquare.edge_size
        + (diagonalSquare.edge_size - 1) * (diagonalSquare.edge_size - 1);
    sort(squares_vec.begin(), squares_vec.end(), compareAreas(diagonalArea));
    return &squares_vec[0];

}



//**********************************************************************************************************
//This part are helper functions for testing and visualization

// print the squares list
void printSquareList(vector<Square>& square_vec) {
    for (auto square_it : square_vec) {
        cout << (square_it);

    }
}
//helper function to check if img and the horizontal square are valids for input img
bool checkSquare(vector<vector<uint8_t>>& img, Square& sq) {
    
    int pos_row = sq.row_corner_point;
    int pos_col = sq.col_corner_point;
    int edge_size = sq.edge_size;
    return img.size()>0 && !img[0].size()>0 &&
        pos_row >= 0 && pos_row < img.size() &&
        pos_col >= 0 && pos_col < img[0].size() &&
        (pos_row + edge_size - 1 < img.size()) && (pos_col + edge_size - 1 < img[0].size());
}
//helper function for building horizontal square in image for testing
void addSquareToImg(vector<vector<uint8_t>>& img, Square& sq) {
    if (!checkSquare(img, sq))
        return;

    int gray_scale = sq.gray_scale_color;
    int pos_row = sq.row_corner_point;
    int pos_col = sq.col_corner_point;
    int edge_size = sq.edge_size;
    // build top_edge
    for (int iter_col = pos_col; iter_col < pos_col + edge_size; ++iter_col) {

        img[pos_row][iter_col] = gray_scale;

    }
    // build bottom_edge
    for (int iter_col = pos_col; iter_col < pos_col + edge_size; ++iter_col) {

        img[pos_row + edge_size - 1][iter_col] = gray_scale;

    }
    // build left_edge
    for (int iter_row = pos_row + 1; iter_row < pos_row + edge_size - 1; ++iter_row) {

        img[iter_row][pos_col] = gray_scale;

    }

    // build right_edge
    for (int iter_row = pos_row + 1; iter_row < pos_row + edge_size - 1; ++iter_row) {

        img[iter_row][pos_col + edge_size - 1] = gray_scale;

    }

}
//helper function to check if img and the diagonal square are valids for input img
bool checkDiagonalSquare(vector<vector<uint8_t>>& img, Square& sq) {

    int pos_row = sq.row_corner_point;
    int pos_col = sq.col_corner_point;
    int edge_size = sq.edge_size;
    return img.size()>0 && img[0].size()>0 &&
        pos_row >= 0 && pos_row < img.size() &&
        pos_col >= 0 && pos_col < img[0].size() &&
        (pos_row + edge_size - 1 < img.size()) && (pos_col + edge_size - 1 < img[0].size()) &&
        (pos_col - (edge_size - 1) < img[0].size()) && (pos_row + 2*(edge_size - 1) < img.size());
}
//helper function for building diagonal square in image for testing
void addDiagonalSquareToImg(vector<vector<uint8_t>>& img, Square& sq) {
    if (!checkDiagonalSquare(img, sq))
        return;
    uint8_t gray_scale = sq.gray_scale_color;
    uint8_t pos_row = sq.row_corner_point;
    uint8_t pos_col = sq.col_corner_point;
    uint8_t edge_size = sq.edge_size;
    uint8_t iter_row;
    uint8_t iter_col;
    // build topRight_edge
    for (iter_row = pos_row, iter_col = pos_col; iter_row < pos_row + edge_size && iter_col < pos_col + edge_size; ++iter_row, ++iter_col) {
        img[iter_row][iter_col] = gray_scale;
    }
    // build topLeft_edge
    for (iter_row = pos_row + 1, iter_col = pos_col - 1; iter_row < pos_row + edge_size && iter_col > pos_col - edge_size; ++iter_row, --iter_col) {
        img[iter_row][iter_col] = gray_scale;
    }
    int new_pos_row = pos_row + edge_size - 1;
    int new_pos_col = pos_col + edge_size - 1;
    // build bottomRight_edge
    for (iter_row = new_pos_row + 1, iter_col = new_pos_col - 1; iter_row < new_pos_row + edge_size && iter_col > new_pos_col - edge_size; ++iter_row, --iter_col) {
        img[iter_row][iter_col] = gray_scale;
    }
    new_pos_col = pos_col - (edge_size - 1);
    // build right_edge
    for (iter_row = new_pos_row + 1, iter_col = new_pos_col + 1; iter_row < new_pos_row + edge_size - 1 && iter_col < new_pos_col + edge_size - 1; ++iter_row, ++iter_col) {
        img[iter_row][iter_col] = gray_scale;
    }
}
// helper function for printing image as a matrix
void printImgMatrix(vector<vector<uint8_t>>& img) {
    for (int row = 0; row < img.size(); ++row) {
        for (int col = 0; col < img[0].size(); ++col) {
            cout << right << setw(4) << +img[row][col] << " ";
        }
        cout << endl;

    }
    cout << endl;
}
//display vector<vector<uint8t>> as an image with openCV
// Note: The display is increase the visability for small pixel size images for better viewing.
void dispalyAsImage(vector<std::vector<uint8_t>>& vec) {
    Mat matDisplay(vec.size(), vec.at(0).size(), CV_8U);
    for (int i = 0; i < matDisplay.rows; ++i)
        for (int j = 0; j < matDisplay.cols; ++j)
            matDisplay.at<uchar>(i, j) = vec.at(i).at(j);
    resize(matDisplay, matDisplay, Size(256, 256), 0, 0,
        INTER_LINEAR);
    imshow("Displaying, close the window to continue or wait 2 sec", matDisplay);
    waitKey(2500);
    /*
    int k = waitKey(0); // for saving image propose
    if (k == 's')
    {
        imwrite("diagonalSquares.png", matDisplay);
    }
    */
    
   
    
}
// visualizing the overlapping of 2 squares (horizontal and diagonal)
// put one square on the other and display
//display the diagonal allways above the horizontal
void visualizeOverLapping(Square& diagonalSquare, Square& horizontalSquare) {
    int diagonalBoardSize = diagonalSquare.edge_size + diagonalSquare.edge_size - 1;
    int horizontalSquareBoardSize = horizontalSquare.edge_size;
    if (diagonalBoardSize > horizontalSquareBoardSize) {
        // most of horizontal square is inside the diagonal square
        vector<uint8_t> zeroes_vec(diagonalBoardSize, 0);

        vector<vector<uint8_t>> board(diagonalBoardSize, zeroes_vec); // building the board

        uint8_t row_diag = 0; //adding the diagonal square
        uint8_t col_diag = diagonalBoardSize / 2;
        Square diag_vis = Square(row_diag, col_diag, diagonalSquare.edge_size, diagonalSquare.gray_scale_color);
        uint8_t row_horiz = row_diag + diagonalSquare.edge_size - horizontalSquare.edge_size / 2; // adding the horizontal square
        uint8_t col_horiz = col_diag - horizontalSquare.edge_size / 2;
        Square horiz_vis = Square(row_horiz, col_horiz, horizontalSquare.edge_size, horizontalSquare.gray_scale_color);
        addSquareToImg(board, horiz_vis);// most of diagonal square is inside the best square
        addDiagonalSquareToImg(board, diag_vis);
        //printImgMatrix(board);
        dispalyAsImage(board);
    }
    else {
        // most of diagonal square is inside the horizontal square
        vector<uint8_t> zeroes_vec(horizontalSquareBoardSize, 0);

        vector<vector<uint8_t>> board(horizontalSquareBoardSize, zeroes_vec); // building the board

        uint8_t row_horiz = 0;//adding the horizontal square
        uint8_t col_horiz = 0;
        Square horiz_vis = Square(row_horiz, col_horiz, horizontalSquare.edge_size, horizontalSquare.gray_scale_color);

        uint8_t row_diag = 0; //adding the diagonal square
        uint8_t col_diag = horizontalSquareBoardSize / 2;
        Square diag_vis = Square(row_diag, col_diag, diagonalSquare.edge_size, diagonalSquare.gray_scale_color);
        addSquareToImg(board, horiz_vis);
        addDiagonalSquareToImg(board, diag_vis);
        //printImgMatrix(board);
        dispalyAsImage(board);
    }
}


//get an image path and return vector
// tranform CV::Mat (from the image) to vector<vector<uint8t>> (using openCV)
vector<std::vector<uchar>> ImageToVector(string path) {
    
	Mat img = imread(path, IMREAD_GRAYSCALE);
	vector<uchar> coloumns_vec(img.cols, 0);
	vector<vector<uchar>> transformed(img.rows, coloumns_vec);
    for (int i = 0; i < img.rows; ++i) {
        for (int j = 0; j < img.cols; ++j) {

            uchar try_val = img.at<uchar>(i, j);
            transformed[i][j] = (uchar)try_val;
        }
    }
	return transformed;
}


//some test in main
int main()
{

    /* if we are getting an image we can transform it to vector (for non open-cv environment)
    using ImageToVector(string path)
    building first image
    row(right top corner) ,col(right top corner)  ,edge_size ,gray_scale_value
    vector<uint8_t> horiz_vec(256, 0);
    vector<vector<uint8_t>> img1(256, horiz_vec);
    */
    /*
    Square s1(4, 150, 50, 250);
    Square s2(50, 50, 20, 250);
    Square s3(25, 25, 20, 249);
    Square s4(170, 170, 80, 248);
    Square s5(240, 0, 10, 247);
    Square s6(0, 20, 15, 246);
    Square s7(200, 100, 30, 246);
    Square s8(90, 30, 30, 246);
    addSquareToImg(img1, s1);
    addSquareToImg(img1, s2);
    addSquareToImg(img1, s3);
    addSquareToImg(img1, s4);
    addSquareToImg(img1, s5);
    addSquareToImg(img1, s6);
    addSquareToImg(img1, s7);
    addSquareToImg(img1, s8);
    //should not add to img
    Square s_lie(10, 10, 250, 70);
    addSquareToImg(img1, s_lie);
    printImgMatrix(im2);
    */
    

    //We built this image befor and now we getting it

    // I can get this 'ImageToVector' function inside 'findAllSquares' function (or continue working with Mat type instead of vector),
    //but for non-opencv enviroment I separtated between thm both
    vector<vector<uchar>> img1 = ImageToVector("horizontalSquares.png");   
    dispalyAsImage(img1);
    vector<Square> ans = findAllSquares(img1);
    cout << "Section 1 output:" << endl;
    printSquareList(ans);
    cout << endl;

    //building second image for diagonal square
    //row(top corner - diagonal) ,col(top corner - diagonal)  ,edge_size ,gray_scale_value
    /*
    vector<uint8_t> horiz_vec2(128, 0);
    vector<vector<uint8_t>> im2_diagonal(128, horiz_vec2);
    Square s_diagonal(10, 60, 10, 245);
    addDiagonalSquareToImg(im2_diagonal, s_diagonal);
    dispalyAsImage(im2_diagonal);
    //printImgMatrix(im2_diagonal);
    */
    
    // I can get this 'ImageToVector' function inside 'findAllSquares' function (or continue working with Mat type instead of vector),
    //but for non-opencv enviroment I separtated between thm both
    vector<vector<uchar>> im2_diagonal = ImageToVector("diagonalSquares.png");
    dispalyAsImage(im2_diagonal);
    //find the diagonal square in the image
    Square* diag_square_ptr = nullptr;
    try {
        diag_square_ptr = findDiagonalSquare(im2_diagonal);

    }
    catch (noSquareException& e) {// no square found
        cout << e.what() << std::endl;
        return 0;
    }

    cout << "Diagonal square properties:" << endl;
    cout << *diag_square_ptr << endl;
    // find the best overlapping square from img1 to diagonal square
    Square* best_fit_horiz_square_ptr = nullptr;
    try {
         best_fit_horiz_square_ptr = findBestOverlap(*diag_square_ptr, ans);

    }
    catch (noSquareException& e)  {// no square found
        cout << e.what() << std::endl;
        return 0;
    }
   
     cout << "Section 3 output (best overfitting square):" << endl;
     cout << *(best_fit_horiz_square_ptr) << endl;
     // display image for visualizing
     cout << "Section 3 visualization (best overfitting square):" << endl;
     visualizeOverLapping(*diag_square_ptr, *best_fit_horiz_square_ptr);
     delete diag_square_ptr;
   
    
	return 0;
}

