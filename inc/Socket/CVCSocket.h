#include <arpa/inet.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>

//*****************
// flags
#define CLEAR 0
#define LOADVOL 1
#define VIEW 2
#define PROJ 4
#define COLTBL 8
#define SHDREN 16
#define TIMEST 32
#define TIMEED 64
#define RESOL 128
#define UPDATEGL 256
#define RENMODE 512
//*****************

#define SERVER 1
#define CLIENT 2
#define BACKLOG 10

#define MAX_DATA_SIZE 1024

//#define DEBUG_SOCKET

typedef struct _msg {
   _msg() {
      bzero(data, MAX_DATA_SIZE );
   }
   short int flag;
   char data[MAX_DATA_SIZE];
   short int size;
}MSG;

class CVCSocketTCP {
public:
   CVCSocketTCP( int flag, int val, int _maxDataSizeByte );
   ~CVCSocketTCP();
   
   bool _connect(void);
   bool _receive(short int *_flag, void* _data, int *_dataSize);
   bool _wreceive(short int *_flag, void* _data, int *_dataSize);
   bool _wreceive( std::vector<MSG> *_msg );
   bool _send( short int _flag, const void* _data, int _dataSize);
   bool setTargets( char *_fname );
   void setTargets( int _nHosts, std::string *_hostnames, int *_ports );

private:
   int m_bufferSize;
   char *m_buffer;
   int m_socket_id;   // for server
   int m_port;
   bool m_connected;
   int *m_socket_idv; // for cleint
   int m_nServer;
   bool *m_connectedv;

   struct sockaddr_in m_clientAddr;
   socklen_t m_clientLen;
   std::vector<struct sockaddr_in> m_targetAddr; // if server, address of nodes in sub call tree
                                               // if client, address of server nodes

   int m_maxDataSizeInByte;

   bool m_socketBinded;
   bool m_hasSubtree;
   int m_IAM;
};

class CVCSocketUDP {
public:
   CVCSocketUDP( int _port, int _maxDataSizeByte );
   ~CVCSocketUDP();

   bool receive(short int *_flag, void* _data, int *_dataSize);
   bool wreceive( short int *_flag, void* _data, int *_dataSize); //wait to receive

   bool send( short int _flag, const void* _data, int _dataSize);
   bool setTargets( char *_fname );
   void setTargets( int _nHosts, std::string *_hostnames, int *_ports );

private:
   int m_bufferSize;
   char *m_buffer;
   int m_socket_id;
   std::vector<struct sockaddr_in> m_targetAddr; // if server, address of nodes in sub call tree
                                               // if client, address of server nodes

   int m_port;
   int m_maxDataSizeInByte;

   bool m_socketBinded;
   bool m_hasSubtree;
};
