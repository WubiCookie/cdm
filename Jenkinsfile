def sendSlack_notification(failed_stage)
{
  def name=powershell(returnStdout: true, script: "git log -1 --pretty=format:'%an'") 
  def url= env.BUILD_URL.replace("localhost","10.37.0.121")
  slackSend (color: "danger", channel: '#ci', message: "Build failed on stage ${failed_stage}: author: ${name} repository: ${env.GIT_URL} commit: ${env.GIT_COMMIT} see error at: ${url}console")
}
pipeline{
  agent any
  triggers {
    pollSCM('0 23 * * *')
  }
  stages{
  stage('build'){
      steps
      {
        script{
          try{
            powershell "xmake -y"         
          }catch(e){
            sendSlack_notification(env.STAGE_NAME)
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
            sendSlack_notification(env.STAGE_NAME)
            throw e
          }
        }  
      }
    }
  }
}
