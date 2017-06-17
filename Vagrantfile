# -*- mode: ruby -*-
# vi: set ft=ruby :

VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
    config.vm.box = "debian/jessie64"

    config.vm.provider "virtualbox" do |vb|
        vb.cpus = 2
        vb.memory = 4096
    end

    config.vm.network :private_network, ip: "192.168.33.12"
    config.vm.network :forwarded_port, host: 8080, guest: 80

    config.vm.synced_folder ".", "/home/vagrant/network_simulation", :nfs => true

    config.vm.provision "shell", inline: "echo 'start provision'"
    config.vm.provision "shell", path: "setup.sh"
    config.vm.provision "shell", inline: "echo 'finish provision'"
end
