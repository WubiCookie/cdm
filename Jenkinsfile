pipeline{
  agent any
  stages{
  stage('build'){
      steps
      {
        script{
          try{
            powershell "xmake -y"         
          }catch(e){
            def name=powershell(returnStdout: true, script: "git log -1 --pretty=format:'%an'") 
            def job_name = env.JOB_NAME
            def table = job_name.split("/")
            def folder=table[0]
            slackSend channel: '#ci', message: "Build failed on stage build: Author: ${name} repository: ${env.GIT_URL} commit: ${env.GIT_COMMIT} see error at: http://10.37.0.121:2020/job/${folder}/job/${env.GIT_BRANCH}/${env.BUILD_ID}/console"
            throw e
          }
        }  
      }
    }
  stage('test'){
      steps
      {
        script{
          try{
            powershell "xmake run"
          }catch(e){
            def name=powershell(returnStdout: true, script: "git log -1 --pretty=format:'%an'") 
            def job_name = env.JOB_NAME
            def table = job_name.split("/")
            def folder=table[0]
            slackSend channel: '#ci', message: "Build failed on stage test: Author: ${name} repository: ${env.GIT_URL} commit: ${env.GIT_COMMIT} see error at: http://10.37.0.121:2020/job/${folder}/job/${env.GIT_BRANCH}/${env.BUILD_ID}/console"
            throw e
          }
        }  
      }
    }
  }
}
